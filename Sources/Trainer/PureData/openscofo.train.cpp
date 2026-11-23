/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include <Trainer.hpp>

#include <atomic>
#include <exception>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <m_pd.h>
#include <g_canvas.h>

enum class OpenScofoTrainMessageType {
    Progress,
    Status,
    Segments,
    Accuracy,
    BalancedAccuracy,
    ClassAccuracy,
    Error,
    Log
};
enum class OpenScofoTrainJob { Analyze, Train, Export, Cut, ClearData };

// ─────────────────────────────────────
struct OpenScofoTrainMessage {
    OpenScofoTrainMessageType Type = OpenScofoTrainMessageType::Status;
    std::string Text;
    double Value = 0.0;
    spdlog::level::level_enum Level = spdlog::level::info;
};

class OpenScofoTrainState;
void openscofo_train_log_callback(const spdlog::details::log_msg &log, void *data);

// ─────────────────────────────────────
class OpenScofoTrainState {
  public:
    OpenScofo::Trainer Trainer;
    OpenScofo::CatConfiguration Config;
    std::mutex QueueMutex;
    std::vector<OpenScofoTrainMessage> Queue;
    std::thread Worker;
    std::atomic<bool> Running{false};
    std::string WorkerStatus;

    explicit OpenScofoTrainState(const std::filesystem::path &patchDirectory) {
        Config.BasePath = patchDirectory;
        Trainer.SetConfiguration(Config);
        Trainer.SetProgressCallback([this](const std::string &status, double progress) {
            if (status != WorkerStatus) {
                WorkerStatus = status;
                Push({OpenScofoTrainMessageType::Status, status, 0.0});
            }
            Push({OpenScofoTrainMessageType::Progress, {}, progress});
        });
        Trainer.SetLogCallback(openscofo_train_log_callback, this);
    }

    ~OpenScofoTrainState() {
        Trainer.Cancel();
        if (Worker.joinable()) {
            Worker.join();
        }
    }

    void Push(OpenScofoTrainMessage message) {
        std::lock_guard<std::mutex> lock(QueueMutex);
        Queue.push_back(std::move(message));
    }

    std::vector<OpenScofoTrainMessage> TakeMessages() {
        std::vector<OpenScofoTrainMessage> messages;
        std::lock_guard<std::mutex> lock(QueueMutex);
        messages.swap(Queue);
        return messages;
    }
};

// ─────────────────────────────────────
void openscofo_train_log_callback(const spdlog::details::log_msg &log, void *data) {
    auto *state = static_cast<OpenScofoTrainState *>(data);
    if (state == nullptr) {
        return;
    }
    state->Push({OpenScofoTrainMessageType::Log, std::string(log.payload.data(), log.payload.size()), 0.0, log.level});
}

t_class *openscofo_train_class = nullptr;

// ─────────────────────────────────────
struct OpenScofoTrain {
    t_object Object;
    OpenScofoTrainState *State;
    t_clock *Clock;
    t_outlet *ProgressOutlet;
    t_outlet *StatusOutlet;
    t_outlet *EvaluationOutlet;
};

// ─────────────────────────────────────
void openscofo_train_output(t_outlet *outlet, const char *selector, const std::string &value) {
    t_atom atom;
    SETSYMBOL(&atom, gensym(value.c_str()));
    outlet_anything(outlet, gensym(selector), 1, &atom);
}

// ─────────────────────────────────────
void openscofo_train_output(t_outlet *outlet, const char *selector, double value) {
    t_atom atom;
    SETFLOAT(&atom, static_cast<t_float>(value));
    outlet_anything(outlet, gensym(selector), 1, &atom);
}

// ─────────────────────────────────────
void openscofo_train_tick(OpenScofoTrain *object) {
    OpenScofoTrainState *state = object->State;
    for (const OpenScofoTrainMessage &message : state->TakeMessages()) {
        switch (message.Type) {
        case OpenScofoTrainMessageType::Progress:
            openscofo_train_output(object->ProgressOutlet, "progress", message.Value);
            break;
        case OpenScofoTrainMessageType::Status:
            openscofo_train_output(object->StatusOutlet, "status", message.Text);
            break;
        case OpenScofoTrainMessageType::Segments:
            openscofo_train_output(object->EvaluationOutlet, "segments", message.Value);
            break;
        case OpenScofoTrainMessageType::Accuracy:
            openscofo_train_output(object->EvaluationOutlet, "accuracy", message.Value);
            break;
        case OpenScofoTrainMessageType::BalancedAccuracy:
            openscofo_train_output(object->EvaluationOutlet, "balanced_accuracy", message.Value);
            break;
        case OpenScofoTrainMessageType::ClassAccuracy: {
            t_atom atoms[2];
            SETSYMBOL(&atoms[0], gensym(message.Text.c_str()));
            SETFLOAT(&atoms[1], static_cast<t_float>(message.Value));
            outlet_anything(object->EvaluationOutlet, gensym("class_accuracy"), 2, atoms);
            break;
        }
        case OpenScofoTrainMessageType::Error:
            openscofo_train_output(object->StatusOutlet, "error", message.Text);
            break;
        case OpenScofoTrainMessageType::Log:
            if (message.Level == spdlog::level::err || message.Level == spdlog::level::critical) {
                logpost(object, PD_ERROR, "[openscofo.train] %s", message.Text.c_str());
            } else if (message.Level == spdlog::level::warn || message.Level == spdlog::level::info) {
                logpost(object, PD_NORMAL, "[openscofo.train] %s", message.Text.c_str());
            } else {
                logpost(object, PD_DEBUG, "[openscofo.train] %s", message.Text.c_str());
            }
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
void openscofo_train_startjob(OpenScofoTrain *object, OpenScofoTrainJob job, std::filesystem::path path = {},
                              std::vector<float> samples = {}, int sampleRate = 0) {
    OpenScofoTrainState *state = object->State;
    if (state->Running.load()) {
        openscofo_train_output(object->StatusOutlet, "error", "openscofo.train is already busy");
        return;
    }
    if (state->Worker.joinable()) {
        state->Worker.join();
    }

    const char *initialStatus = "analyzing";
    if (job == OpenScofoTrainJob::Train)
        initialStatus = "training";
    else if (job == OpenScofoTrainJob::Export)
        initialStatus = "exporting";
    else if (job == OpenScofoTrainJob::Cut)
        initialStatus = "cutting";
    else if (job == OpenScofoTrainJob::ClearData)
        initialStatus = "clearing";

    state->WorkerStatus = initialStatus;
    state->Running.store(true);
    openscofo_train_output(object->StatusOutlet, "status", initialStatus);
    state->Worker =
        std::thread([state, job, path = std::move(path), samples = std::move(samples), sampleRate]() mutable {
            bool succeeded = false;
            try {
                if (job == OpenScofoTrainJob::Analyze) {
                    succeeded = state->Trainer.Analyze();
                } else if (job == OpenScofoTrainJob::Train) {
                    succeeded = state->Trainer.Train();
                    const OpenScofo::Evaluation &evaluation = state->Trainer.GetEvaluation();
                    if (succeeded && evaluation.TestedSamples != 0) {
                        state->Push({OpenScofoTrainMessageType::Accuracy, {}, evaluation.Accuracy});
                        state->Push({OpenScofoTrainMessageType::BalancedAccuracy, {}, evaluation.BalancedAccuracy});
                        for (const auto &[label, accuracy] : evaluation.ClassAccuracy) {
                            state->Push({OpenScofoTrainMessageType::ClassAccuracy, label, accuracy});
                        }
                    }
                } else if (job == OpenScofoTrainJob::Export) {
                    succeeded = state->Trainer.ExportModel(path);
                } else if (job == OpenScofoTrainJob::Cut) {
                    const std::optional<std::size_t> exported = state->Trainer.Cut(samples, sampleRate, path);
                    succeeded = exported.has_value();
                    if (exported) {
                        state->Push({OpenScofoTrainMessageType::Segments, {}, static_cast<double>(*exported)});
                    }
                } else {
                    state->Trainer.ClearData();
                    succeeded = true;
                }
            } catch (const std::exception &error) {
                state->Push({OpenScofoTrainMessageType::Error, error.what(), 0.0});
            }
            if (succeeded) {
                state->Push({OpenScofoTrainMessageType::Status, "done", 0.0});
            } else if (state->Trainer.IsCancelled()) {
                state->Push({OpenScofoTrainMessageType::Status, "cancelled", 0.0});
            } else {
                state->Push({OpenScofoTrainMessageType::Status, "failed", 0.0});
            }
            state->Running.store(false);
        });
    clock_delay(object->Clock, 0.0);
}

// ─────────────────────────────────────
bool openscofo_train_canconfigure(OpenScofoTrain *object) {
    if (object->State->Running.load()) {
        openscofo_train_output(object->StatusOutlet, "error", "cannot change settings while a job is running");
        return false;
    }
    if (object->State->Worker.joinable()) {
        object->State->Worker.join();
    }
    return true;
}

// ─────────────────────────────────────
void openscofo_train_applyconfiguration(OpenScofoTrain *object) {
    object->State->Trainer.SetConfiguration(object->State->Config);
}

// ─────────────────────────────────────
bool openscofo_train_requiresymbol(OpenScofoTrain *object, int argc, t_atom *argv, const char *message) {
    if (argc != 1 || argv[0].a_type != A_SYMBOL) {
        openscofo_train_output(object->StatusOutlet, "error", message);
        return false;
    }
    return true;
}

// ─────────────────────────────────────
void openscofo_train_setfolder(OpenScofoTrain *object, t_symbol *, int argc, t_atom *argv) {
    if (!openscofo_train_canconfigure(object) ||
        !openscofo_train_requiresymbol(object, argc, argv, "folder requires one path"))
        return;
    object->State->Config.DatasetFolder = atom_getsymbol(argv)->s_name;
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setirfolders(OpenScofoTrain *object, t_symbol *, int argc, t_atom *argv) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.ImpulseResponseFolders.clear();
    for (int i = 0; i < argc; ++i) {
        if (argv[i].a_type != A_SYMBOL) {
            openscofo_train_output(object->StatusOutlet, "error", "ir_folder paths must be symbols");
            return;
        }
        object->State->Config.ImpulseResponseFolders.emplace_back(atom_getsymbol(argv + i)->s_name);
    }
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setdescriptors(OpenScofoTrain *object, t_symbol *, int argc, t_atom *argv) {
    if (!openscofo_train_canconfigure(object))
        return;
    if (argc == 0) {
        openscofo_train_output(object->StatusOutlet, "error", "descriptors requires at least one name");
        return;
    }
    object->State->Config.Descriptors.clear();
    for (int i = 0; i < argc; ++i) {
        if (argv[i].a_type != A_SYMBOL) {
            openscofo_train_output(object->StatusOutlet, "error", "descriptor names must be symbols");
            return;
        }
        object->State->Config.Descriptors.emplace_back(atom_getsymbol(argv + i)->s_name);
    }
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_settrainingsamplerates(OpenScofoTrain *object, t_symbol *, int argc, t_atom *argv) {
    if (!openscofo_train_canconfigure(object))
        return;
    if (argc == 0) {
        openscofo_train_output(object->StatusOutlet, "error", "training_sample_rates requires at least one rate");
        return;
    }
    object->State->Config.TrainingSampleRates.clear();
    for (int i = 0; i < argc; ++i) {
        if (argv[i].a_type != A_FLOAT) {
            openscofo_train_output(object->StatusOutlet, "error", "training sample rates must be numbers");
            return;
        }
        object->State->Config.TrainingSampleRates.push_back(atom_getint(argv + i));
    }
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setsamplerate(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.SampleRate = static_cast<int>(value);
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setfftsize(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.FFTSize = static_cast<int>(value);
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_sethopsize(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.HopSize = static_cast<int>(value);
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setmindb(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.MinDB = value;
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setmaxpitchconfidence(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.MaxPitchConfidence = value;
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setaugmentationfactor(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.AugmentationFactor = value;
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_settestfraction(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.TestFraction = value;
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setiterations(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.Iterations = static_cast<int>(value);
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setlearningrate(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.LearningRate = value;
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setdepth(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.Depth = static_cast<int>(value);
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setearlystopping(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.EarlyStoppingRounds = static_cast<int>(value);
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_setseed(OpenScofoTrain *object, t_float value) {
    if (!openscofo_train_canconfigure(object))
        return;
    object->State->Config.RandomSeed = static_cast<int>(value);
    openscofo_train_applyconfiguration(object);
}

// ─────────────────────────────────────
void openscofo_train_analyze(OpenScofoTrain *object) {
    openscofo_train_startjob(object, OpenScofoTrainJob::Analyze);
}

// ─────────────────────────────────────
void openscofo_train_train(OpenScofoTrain *object) {
    openscofo_train_startjob(object, OpenScofoTrainJob::Train);
}

// ─────────────────────────────────────
void openscofo_train_export(OpenScofoTrain *object, t_symbol *path) {
    if (path == nullptr || path == &s_ || path->s_name[0] == '\0') {
        openscofo_train_output(object->StatusOutlet, "error", "export requires one model path");
        return;
    }
    openscofo_train_startjob(object, OpenScofoTrainJob::Export, path->s_name);
}

// ─────────────────────────────────────
void openscofo_train_cut(OpenScofoTrain *object, t_symbol *, int argc, t_atom *argv) {
    if (!openscofo_train_canconfigure(object)) {
        return;
    }
    if (argc != 2 || argv[0].a_type != A_SYMBOL || argv[1].a_type != A_SYMBOL) {
        openscofo_train_output(object->StatusOutlet, "error", "cut requires an array and a label");
        return;
    }

    t_symbol *arrayName = atom_getsymbol(argv);
    t_garray *array = reinterpret_cast<t_garray *>(pd_findbyclass(arrayName, garray_class));
    if (array == nullptr) {
        openscofo_train_output(object->StatusOutlet, "error", "array was not found");
        return;
    }
    int size = 0;
    t_word *words = nullptr;
    if (!garray_getfloatwords(array, &size, &words) || size <= 0 || words == nullptr) {
        openscofo_train_output(object->StatusOutlet, "error", "array does not contain float samples");
        return;
    }

    std::vector<float> samples(static_cast<std::size_t>(size));
    for (int index = 0; index < size; ++index) {
        samples[static_cast<std::size_t>(index)] = words[index].w_float;
    }
    const int sampleRate = sys_getsr() > 0 ? static_cast<int>(sys_getsr()) : object->State->Config.SampleRate;
    openscofo_train_startjob(object, OpenScofoTrainJob::Cut, atom_getsymbol(argv + 1)->s_name, std::move(samples),
                             sampleRate);
}

// ─────────────────────────────────────
void openscofo_train_cleardata(OpenScofoTrain *object) {
    openscofo_train_startjob(object, OpenScofoTrainJob::ClearData);
}

// ─────────────────────────────────────
void openscofo_train_cancel(OpenScofoTrain *object) {
    if (!object->State->Running.load()) {
        openscofo_train_output(object->StatusOutlet, "status", "idle");
        return;
    }
    object->State->Trainer.Cancel();
    openscofo_train_output(object->StatusOutlet, "status", "cancelling");
}

// ─────────────────────────────────────
void *openscofo_train_new(t_symbol *, int, t_atom *) {
    auto *object = reinterpret_cast<OpenScofoTrain *>(pd_new(openscofo_train_class));
    if (object == nullptr) {
        return nullptr;
    }
    t_canvas *canvas = canvas_getcurrent();
    const std::filesystem::path patchDirectory = canvas_getdir(canvas)->s_name;
    try {
        object->State = new OpenScofoTrainState(patchDirectory);
    } catch (const std::exception &error) {
        pd_error(object, "[openscofo.train] %s", error.what());
        return nullptr;
    }
    object->ProgressOutlet = outlet_new(&object->Object, &s_list);
    object->StatusOutlet = outlet_new(&object->Object, &s_list);
    object->EvaluationOutlet = outlet_new(&object->Object, &s_list);
    object->Clock = clock_new(object, reinterpret_cast<t_method>(openscofo_train_tick));
    return object;
}

// ─────────────────────────────────────
void openscofo_train_free(OpenScofoTrain *object) {
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

// ─────────────────────────────────────
#ifdef _WIN32
extern "C" __declspec(dllexport) void setup_openscofo0x2etrain(void) {
#else
extern "C" __attribute__((visibility("default"))) void setup_openscofo0x2etrain(void) {
#endif
    openscofo_train_class =
        class_new(gensym("openscofo.train"), reinterpret_cast<t_newmethod>(openscofo_train_new),
                  reinterpret_cast<t_method>(openscofo_train_free), sizeof(OpenScofoTrain), CLASS_DEFAULT, A_GIMME, 0);

    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setfolder), gensym("folder"),
                    A_GIMME, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setirfolders),
                    gensym("ir_folder"), A_GIMME, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setdescriptors),
                    gensym("descriptors"), A_GIMME, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setsamplerate),
                    gensym("sample_rate"), A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_settrainingsamplerates),
                    gensym("training_sample_rates"), A_GIMME, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setfftsize), gensym("fft_size"),
                    A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_sethopsize), gensym("hop_size"),
                    A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setmindb), gensym("min_db"),
                    A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setmaxpitchconfidence),
                    gensym("max_pitch_confidence"), A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setaugmentationfactor),
                    gensym("augmentation_factor"), A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_settestfraction),
                    gensym("test_fraction"), A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setiterations),
                    gensym("iterations"), A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setlearningrate),
                    gensym("learning_rate"), A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setdepth), gensym("depth"),
                    A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setearlystopping),
                    gensym("early_stopping_rounds"), A_FLOAT, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_setseed), gensym("seed"), A_FLOAT,
                    0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_analyze), gensym("analyze"),
                    A_NULL);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_train), gensym("train"), A_NULL);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_export), gensym("export"),
                    A_SYMBOL, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_cut), gensym("cut"), A_GIMME, 0);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_cleardata), gensym("clear_data"),
                    A_NULL);
    class_addmethod(openscofo_train_class, reinterpret_cast<t_method>(openscofo_train_cancel), gensym("cancel"),
                    A_NULL);

    logpost(nullptr, PD_NORMAL, "[openscofo.train] native trainer %s", OPENSCOFO_VERSION);
}
