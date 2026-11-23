/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include <Trainer.hpp>

#include <atomic>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <m_pd.h>
#include <g_canvas.h>

namespace {

enum class MessageType { Progress, Status, Accuracy, BalancedAccuracy, ClassAccuracy, Error, Log };
enum class Job { Analyze, Train, Export, ClearData };

// ─────────────────────────────────────
struct Message {
    MessageType Type = MessageType::Status;
    std::string Text;
    double Value = 0.0;
};

// ─────────────────────────────────────
class PdTrainerState {
  public:
    OpenScofo::Trainer Trainer;
    OpenScofo::CatConfiguration Config;
    std::mutex QueueMutex;
    std::vector<Message> Queue;
    std::thread Worker;
    std::atomic<bool> Running{false};
    std::string WorkerStatus;

    explicit PdTrainerState(const std::filesystem::path &patchDirectory) {
        Config.BasePath = patchDirectory;
        Trainer.SetConfiguration(Config);
        Trainer.SetProgressCallback([this](const std::string &status, double progress) {
            if (status != WorkerStatus) {
                WorkerStatus = status;
                Push({MessageType::Status, status, 0.0});
            }
            Push({MessageType::Progress, {}, progress});
        });
        // Trainer.SetLogCallback([this](const std::string &message) { Push({MessageType::Log, message, 0.0}); });
    }

    ~PdTrainerState() {
        Trainer.Cancel();
        if (Worker.joinable()) {
            Worker.join();
        }
    }

    void Push(Message message) {
        std::lock_guard<std::mutex> lock(QueueMutex);
        Queue.push_back(std::move(message));
    }

    std::vector<Message> TakeMessages() {
        std::vector<Message> messages;
        std::lock_guard<std::mutex> lock(QueueMutex);
        messages.swap(Queue);
        return messages;
    }
};

t_class *TrainerClass = nullptr;

// ─────────────────────────────────────
struct PdTrainer {
    t_object Object;
    PdTrainerState *State;
    t_clock *Clock;
    t_outlet *ProgressOutlet;
    t_outlet *StatusOutlet;
    t_outlet *EvaluationOutlet;
};

// ─────────────────────────────────────
void OutputSymbol(t_outlet *outlet, const char *selector, const std::string &value) {
    t_atom atom;
    SETSYMBOL(&atom, gensym(value.c_str()));
    outlet_anything(outlet, gensym(selector), 1, &atom);
}

// ─────────────────────────────────────
void OutputFloat(t_outlet *outlet, const char *selector, double value) {
    t_atom atom;
    SETFLOAT(&atom, static_cast<t_float>(value));
    outlet_anything(outlet, gensym(selector), 1, &atom);
}

// ─────────────────────────────────────
void TrainerTick(PdTrainer *object) {
    PdTrainerState *state = object->State;
    for (const Message &message : state->TakeMessages()) {
        switch (message.Type) {
        case MessageType::Progress:
            OutputFloat(object->ProgressOutlet, "progress", message.Value);
            break;
        case MessageType::Status:
            OutputSymbol(object->StatusOutlet, "status", message.Text);
            break;
        case MessageType::Accuracy:
            OutputFloat(object->EvaluationOutlet, "accuracy", message.Value);
            break;
        case MessageType::BalancedAccuracy:
            OutputFloat(object->EvaluationOutlet, "balanced_accuracy", message.Value);
            break;
        case MessageType::ClassAccuracy: {
            t_atom atoms[2];
            SETSYMBOL(&atoms[0], gensym(message.Text.c_str()));
            SETFLOAT(&atoms[1], static_cast<t_float>(message.Value));
            outlet_anything(object->EvaluationOutlet, gensym("class_accuracy"), 2, atoms);
            break;
        }
        case MessageType::Error:
            OutputSymbol(object->StatusOutlet, "error", message.Text);
            break;
        case MessageType::Log:
            logpost(object, PD_NORMAL, "[openscofo.train] %s", message.Text.c_str());
            break;
        }
    }

    if (state->Running.load()) {
        clock_delay(object->Clock, 20.0);
    } else if (state->Worker.joinable()) {
        state->Worker.join();
        // A worker can enqueue its final message after TakeMessages() but
        // before Running becomes false.  One final tick drains that race.
        clock_delay(object->Clock, 0.0);
    }
}

// ─────────────────────────────────────

void StartJob(PdTrainer *object, Job job, std::filesystem::path exportPath = {}) {
    PdTrainerState *state = object->State;
    if (state->Running.load()) {
        OutputSymbol(object->StatusOutlet, "error", "openscofo.train is already busy");
        return;
    }
    if (state->Worker.joinable()) {
        state->Worker.join();
    }

    const char *initialStatus = "analyzing";
    if (job == Job::Train)
        initialStatus = "training";
    else if (job == Job::Export)
        initialStatus = "exporting";
    else if (job == Job::ClearData)
        initialStatus = "clearing";

    state->WorkerStatus = initialStatus;
    state->Running.store(true);
    OutputSymbol(object->StatusOutlet, "status", initialStatus);
    state->Worker = std::thread([state, job, exportPath = std::move(exportPath)]() {
        try {
            if (job == Job::Analyze) {
                state->Trainer.Analyze();
            } else if (job == Job::Train) {
                state->Trainer.Train();
                const OpenScofo::Evaluation &evaluation = state->Trainer.GetEvaluation();
                if (evaluation.TestedSamples != 0) {
                    state->Push({MessageType::Accuracy, {}, evaluation.Accuracy});
                    state->Push({MessageType::BalancedAccuracy, {}, evaluation.BalancedAccuracy});
                    for (const auto &[label, accuracy] : evaluation.ClassAccuracy) {
                        state->Push({MessageType::ClassAccuracy, label, accuracy});
                    }
                }
            } else if (job == Job::Export) {
                state->Trainer.ExportModel(exportPath);
            } else {
                state->Trainer.ClearData();
            }
            state->Push({MessageType::Status, "done", 0.0});
        } catch (const std::exception &error) {
            const std::string message = error.what();
            if (message == "Training cancelled") {
                state->Push({MessageType::Status, "cancelled", 0.0});
            } else {
                state->Push({MessageType::Error, message, 0.0});
            }
        }
        state->Running.store(false);
    });
    clock_delay(object->Clock, 0.0);
}

// ─────────────────────────────────────
bool CanConfigure(PdTrainer *object) {
    if (object->State->Running.load()) {
        OutputSymbol(object->StatusOutlet, "error", "cannot change settings while a job is running");
        return false;
    }
    if (object->State->Worker.joinable()) {
        object->State->Worker.join();
    }
    return true;
}

// ─────────────────────────────────────
void ApplyConfiguration(PdTrainer *object) {
    object->State->Trainer.SetConfiguration(object->State->Config);
}

// ─────────────────────────────────────
bool RequireSymbol(PdTrainer *object, int argc, t_atom *argv, const char *message) {
    if (argc != 1 || argv[0].a_type != A_SYMBOL) {
        OutputSymbol(object->StatusOutlet, "error", message);
        return false;
    }
    return true;
}

// ─────────────────────────────────────
void SetFolder(PdTrainer *object, t_symbol *, int argc, t_atom *argv) {
    if (!CanConfigure(object) || !RequireSymbol(object, argc, argv, "folder requires one path"))
        return;
    object->State->Config.DatasetFolder = atom_getsymbol(argv)->s_name;
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetIRFolders(PdTrainer *object, t_symbol *, int argc, t_atom *argv) {
    if (!CanConfigure(object))
        return;
    object->State->Config.ImpulseResponseFolders.clear();
    for (int i = 0; i < argc; ++i) {
        if (argv[i].a_type != A_SYMBOL) {
            OutputSymbol(object->StatusOutlet, "error", "ir_folder paths must be symbols");
            return;
        }
        object->State->Config.ImpulseResponseFolders.emplace_back(atom_getsymbol(argv + i)->s_name);
    }
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetDescriptors(PdTrainer *object, t_symbol *, int argc, t_atom *argv) {
    if (!CanConfigure(object))
        return;
    if (argc == 0) {
        OutputSymbol(object->StatusOutlet, "error", "descriptors requires at least one name");
        return;
    }
    object->State->Config.Descriptors.clear();
    for (int i = 0; i < argc; ++i) {
        if (argv[i].a_type != A_SYMBOL) {
            OutputSymbol(object->StatusOutlet, "error", "descriptor names must be symbols");
            return;
        }
        object->State->Config.Descriptors.emplace_back(atom_getsymbol(argv + i)->s_name);
    }
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetTrainingSampleRates(PdTrainer *object, t_symbol *, int argc, t_atom *argv) {
    if (!CanConfigure(object))
        return;
    if (argc == 0) {
        OutputSymbol(object->StatusOutlet, "error", "training_sample_rates requires at least one rate");
        return;
    }
    object->State->Config.TrainingSampleRates.clear();
    for (int i = 0; i < argc; ++i) {
        if (argv[i].a_type != A_FLOAT) {
            OutputSymbol(object->StatusOutlet, "error", "training sample rates must be numbers");
            return;
        }
        object->State->Config.TrainingSampleRates.push_back(atom_getint(argv + i));
    }
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetSampleRate(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.SampleRate = static_cast<int>(value);
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetFFTSize(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.FFTSize = static_cast<int>(value);
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetHopSize(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.HopSize = static_cast<int>(value);
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetMinDB(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.MinDB = value;
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetMaxPitchConfidence(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.MaxPitchConfidence = value;
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetAugmentationFactor(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.AugmentationFactor = value;
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetTestFraction(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.TestFraction = value;
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetIterations(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.Iterations = static_cast<int>(value);
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetLearningRate(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.LearningRate = value;
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetDepth(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.Depth = static_cast<int>(value);
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetEarlyStopping(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.EarlyStoppingRounds = static_cast<int>(value);
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void SetSeed(PdTrainer *object, t_float value) {
    if (!CanConfigure(object))
        return;
    object->State->Config.RandomSeed = static_cast<int>(value);
    ApplyConfiguration(object);
}

// ─────────────────────────────────────
void Analyze(PdTrainer *object) {
    StartJob(object, Job::Analyze);
}

// ─────────────────────────────────────
void Train(PdTrainer *object) {
    StartJob(object, Job::Train);
}

// ─────────────────────────────────────
void Export(PdTrainer *object, t_symbol *path) {
    if (path == nullptr || path == &s_ || path->s_name[0] == '\0') {
        OutputSymbol(object->StatusOutlet, "error", "export requires one model path");
        return;
    }
    StartJob(object, Job::Export, path->s_name);
}

// ─────────────────────────────────────
void ClearData(PdTrainer *object) {
    StartJob(object, Job::ClearData);
}

// ─────────────────────────────────────
void Cancel(PdTrainer *object) {
    if (!object->State->Running.load()) {
        OutputSymbol(object->StatusOutlet, "status", "idle");
        return;
    }
    object->State->Trainer.Cancel();
    OutputSymbol(object->StatusOutlet, "status", "cancelling");
}

// ─────────────────────────────────────
void *NewTrainer(t_symbol *, int, t_atom *) {
    PdTrainer *object = reinterpret_cast<PdTrainer *>(pd_new(TrainerClass));
    if (object == nullptr) {
        return nullptr;
    }
    t_canvas *canvas = canvas_getcurrent();
    const std::filesystem::path patchDirectory = canvas_getdir(canvas)->s_name;
    try {
        object->State = new PdTrainerState(patchDirectory);
    } catch (const std::exception &error) {
        pd_error(object, "[openscofo.train] %s", error.what());
        return nullptr;
    }
    object->ProgressOutlet = outlet_new(&object->Object, &s_list);
    object->StatusOutlet = outlet_new(&object->Object, &s_list);
    object->EvaluationOutlet = outlet_new(&object->Object, &s_list);
    object->Clock = clock_new(object, reinterpret_cast<t_method>(TrainerTick));
    return object;
}

// ─────────────────────────────────────
void FreeTrainer(PdTrainer *object) {
    if (object->Clock != nullptr) {
        clock_unset(object->Clock);
    }
    delete object->State;
    object->State = nullptr;
    if (object->Clock != nullptr) {
        clock_free(object->Clock);
        object->Clock = nullptr;
    }
}

} // namespace

// ─────────────────────────────────────
#ifdef _WIN32
extern "C" __declspec(dllexport) void setup_openscofo0x2etrain(void) {
#else
extern "C" __attribute__((visibility("default"))) void setup_openscofo0x2etrain(void) {
#endif
    TrainerClass = class_new(gensym("openscofo.train"), reinterpret_cast<t_newmethod>(NewTrainer),
                             reinterpret_cast<t_method>(FreeTrainer), sizeof(PdTrainer), CLASS_DEFAULT, A_GIMME, 0);

    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetFolder), gensym("folder"), A_GIMME, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetIRFolders), gensym("ir_folder"), A_GIMME, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetDescriptors), gensym("descriptors"), A_GIMME, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetSampleRate), gensym("sample_rate"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetTrainingSampleRates), gensym("training_sample_rates"),
                    A_GIMME, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetFFTSize), gensym("fft_size"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetHopSize), gensym("hop_size"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetMinDB), gensym("min_db"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetMaxPitchConfidence), gensym("max_pitch_confidence"),
                    A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetAugmentationFactor), gensym("augmentation_factor"),
                    A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetTestFraction), gensym("test_fraction"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetIterations), gensym("iterations"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetLearningRate), gensym("learning_rate"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetDepth), gensym("depth"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetEarlyStopping), gensym("early_stopping_rounds"),
                    A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(SetSeed), gensym("seed"), A_FLOAT, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(Analyze), gensym("analyze"), A_NULL);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(Train), gensym("train"), A_NULL);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(Export), gensym("export"), A_SYMBOL, 0);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(ClearData), gensym("clear_data"), A_NULL);
    class_addmethod(TrainerClass, reinterpret_cast<t_method>(Cancel), gensym("cancel"), A_NULL);

    logpost(nullptr, PD_NORMAL, "[openscofo.train] native trainer %s", OPENSCOFO_VERSION);
}
