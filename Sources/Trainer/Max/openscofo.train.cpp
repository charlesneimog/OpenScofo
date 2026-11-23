/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include <Trainer.hpp>

#include <ext.h>
#include <ext_buffer.h>
#include <ext_path.h>

#include <atomic>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <limits>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <utility>
#include <vector>

enum class MaxTrainMessageType { Progress, Status, Segments, Accuracy, BalancedAccuracy, ClassAccuracy, Error, Log };
enum class MaxTrainJob { Analyze, Train, Export, Cut, ClearData };

struct MaxTrainMessage {
    MaxTrainMessageType Type = MaxTrainMessageType::Status;
    std::string Text;
    double Value = 0.0;
    spdlog::level::level_enum Level = spdlog::level::info;
};

class MaxTrainState;
static void max_train_log_callback(const spdlog::details::log_msg &log, void *data);

class MaxTrainState {
  public:
    OpenScofo::Trainer Trainer;
    OpenScofo::CatConfiguration Config;
    std::mutex QueueMutex;
    std::vector<MaxTrainMessage> Queue;
    std::thread Worker;
    std::atomic<bool> Running{false};
    std::string WorkerStatus;

    explicit MaxTrainState(const std::filesystem::path &patchDirectory) {
        Config.BasePath = patchDirectory;
        Trainer.SetConfiguration(Config);
        Trainer.SetProgressCallback([this](const std::string &status, double progress) {
            if (status != WorkerStatus) {
                WorkerStatus = status;
                Push({MaxTrainMessageType::Status, status, 0.0});
            }
            Push({MaxTrainMessageType::Progress, {}, progress});
        });
        Trainer.SetLogCallback(max_train_log_callback, this);
    }

    ~MaxTrainState() {
        Trainer.Cancel();
        if (Worker.joinable()) {
            Worker.join();
        }
    }

    void Push(MaxTrainMessage message) {
        std::lock_guard<std::mutex> lock(QueueMutex);
        Queue.push_back(std::move(message));
    }

    std::vector<MaxTrainMessage> TakeMessages() {
        std::vector<MaxTrainMessage> messages;
        std::lock_guard<std::mutex> lock(QueueMutex);
        messages.swap(Queue);
        return messages;
    }
};

static void max_train_log_callback(const spdlog::details::log_msg &log, void *data) {
    auto *state = static_cast<MaxTrainState *>(data);
    if (state != nullptr) {
        state->Push({MaxTrainMessageType::Log, std::string(log.payload.data(), log.payload.size()), 0.0, log.level});
    }
}

struct MaxOpenScofoTrain {
    t_object Object;
    MaxTrainState *State;
    t_clock *Clock;
    void *ProgressOutlet;
    void *StatusOutlet;
    void *EvaluationOutlet;
};

static t_class *max_train_class = nullptr;

static void max_train_output(void *outlet, const char *selector, const std::string &value) {
    t_atom atom;
    atom_setsym(&atom, gensym(value.c_str()));
    outlet_anything(outlet, gensym(selector), 1, &atom);
}

static void max_train_output(void *outlet, const char *selector, double value) {
    t_atom atom;
    atom_setfloat(&atom, value);
    outlet_anything(outlet, gensym(selector), 1, &atom);
}

static void max_train_tick(MaxOpenScofoTrain *object) {
    MaxTrainState *state = object->State;
    for (const MaxTrainMessage &message : state->TakeMessages()) {
        switch (message.Type) {
        case MaxTrainMessageType::Progress:
            max_train_output(object->ProgressOutlet, "progress", message.Value);
            break;
        case MaxTrainMessageType::Status:
            max_train_output(object->StatusOutlet, "status", message.Text);
            break;
        case MaxTrainMessageType::Segments:
            max_train_output(object->EvaluationOutlet, "segments", message.Value);
            break;
        case MaxTrainMessageType::Accuracy:
            max_train_output(object->EvaluationOutlet, "accuracy", message.Value);
            break;
        case MaxTrainMessageType::BalancedAccuracy:
            max_train_output(object->EvaluationOutlet, "balanced_accuracy", message.Value);
            break;
        case MaxTrainMessageType::ClassAccuracy: {
            t_atom atoms[2];
            atom_setsym(atoms, gensym(message.Text.c_str()));
            atom_setfloat(atoms + 1, message.Value);
            outlet_anything(object->EvaluationOutlet, gensym("class_accuracy"), 2, atoms);
            break;
        }
        case MaxTrainMessageType::Error:
            max_train_output(object->StatusOutlet, "error", message.Text);
            break;
        case MaxTrainMessageType::Log:
            if (message.Level == spdlog::level::err || message.Level == spdlog::level::critical) {
                object_error(reinterpret_cast<t_object *>(object), "%s", message.Text.c_str());
            } else {
                object_post(reinterpret_cast<t_object *>(object), "%s", message.Text.c_str());
            }
            break;
        }
    }

    if (state->Running.load()) {
        clock_delay(object->Clock, 20.0);
    } else if (state->Worker.joinable()) {
        state->Worker.join();
        clock_delay(object->Clock, 0.0);
    }
}

static void max_train_start_job(MaxOpenScofoTrain *object, MaxTrainJob job, std::filesystem::path path = {},
                                std::vector<float> samples = {}, int sampleRate = 0) {
    MaxTrainState *state = object->State;
    if (state->Running.load()) {
        max_train_output(object->StatusOutlet, "error", "openscofo.train is already busy");
        return;
    }
    if (state->Worker.joinable()) {
        state->Worker.join();
    }

    const char *initialStatus = "analyzing";
    if (job == MaxTrainJob::Train)
        initialStatus = "training";
    else if (job == MaxTrainJob::Export)
        initialStatus = "exporting";
    else if (job == MaxTrainJob::Cut)
        initialStatus = "cutting";
    else if (job == MaxTrainJob::ClearData)
        initialStatus = "clearing";

    state->WorkerStatus = initialStatus;
    state->Running.store(true);
    max_train_output(object->StatusOutlet, "status", initialStatus);
    state->Worker =
        std::thread([state, job, path = std::move(path), samples = std::move(samples), sampleRate]() mutable {
            bool succeeded = false;
            try {
                if (job == MaxTrainJob::Analyze) {
                    succeeded = state->Trainer.Analyze();
                } else if (job == MaxTrainJob::Train) {
                    succeeded = state->Trainer.Train();
                    const OpenScofo::Evaluation &evaluation = state->Trainer.GetEvaluation();
                    if (succeeded && evaluation.TestedSamples != 0) {
                        state->Push({MaxTrainMessageType::Accuracy, {}, evaluation.Accuracy});
                        state->Push({MaxTrainMessageType::BalancedAccuracy, {}, evaluation.BalancedAccuracy});
                        for (const auto &[label, accuracy] : evaluation.ClassAccuracy) {
                            state->Push({MaxTrainMessageType::ClassAccuracy, label, accuracy});
                        }
                    }
                } else if (job == MaxTrainJob::Export) {
                    succeeded = state->Trainer.ExportModel(path);
                } else if (job == MaxTrainJob::Cut) {
                    const std::optional<std::size_t> exported = state->Trainer.Cut(samples, sampleRate, path);
                    succeeded = exported.has_value();
                    if (exported) {
                        state->Push({MaxTrainMessageType::Segments, {}, static_cast<double>(*exported)});
                    }
                } else {
                    state->Trainer.ClearData();
                    succeeded = true;
                }
            } catch (const std::exception &error) {
                state->Push({MaxTrainMessageType::Error, error.what(), 0.0});
            }
            if (succeeded) {
                state->Push({MaxTrainMessageType::Status, "done", 0.0});
            } else if (state->Trainer.IsCancelled()) {
                state->Push({MaxTrainMessageType::Status, "cancelled", 0.0});
            } else {
                state->Push({MaxTrainMessageType::Status, "failed", 0.0});
            }
            state->Running.store(false);
        });
    clock_delay(object->Clock, 0.0);
}

static bool max_train_can_configure(MaxOpenScofoTrain *object) {
    if (object->State->Running.load()) {
        max_train_output(object->StatusOutlet, "error", "cannot change settings while a job is running");
        return false;
    }
    if (object->State->Worker.joinable()) {
        object->State->Worker.join();
    }
    return true;
}

static void max_train_apply_configuration(MaxOpenScofoTrain *object) {
    object->State->Trainer.SetConfiguration(object->State->Config);
}

static bool max_train_is_number(const t_atom &atom) {
    return atom.a_type == A_LONG || atom.a_type == A_FLOAT;
}

static void max_train_set_folder(MaxOpenScofoTrain *object, t_symbol *, long argc, t_atom *argv) {
    if (!max_train_can_configure(object))
        return;
    if (argc != 1 || atom_gettype(argv) != A_SYM) {
        max_train_output(object->StatusOutlet, "error", "folder requires one path");
        return;
    }
    object->State->Config.DatasetFolder = atom_getsym(argv)->s_name;
    max_train_apply_configuration(object);
}

static void max_train_set_ir_folders(MaxOpenScofoTrain *object, t_symbol *, long argc, t_atom *argv) {
    if (!max_train_can_configure(object))
        return;
    object->State->Config.ImpulseResponseFolders.clear();
    for (long index = 0; index < argc; ++index) {
        if (atom_gettype(argv + index) != A_SYM) {
            max_train_output(object->StatusOutlet, "error", "ir_folder paths must be symbols");
            return;
        }
        object->State->Config.ImpulseResponseFolders.emplace_back(atom_getsym(argv + index)->s_name);
    }
    max_train_apply_configuration(object);
}

static void max_train_set_descriptors(MaxOpenScofoTrain *object, t_symbol *, long argc, t_atom *argv) {
    if (!max_train_can_configure(object))
        return;
    if (argc == 0) {
        max_train_output(object->StatusOutlet, "error", "descriptors requires at least one name");
        return;
    }
    object->State->Config.Descriptors.clear();
    for (long index = 0; index < argc; ++index) {
        if (atom_gettype(argv + index) != A_SYM) {
            max_train_output(object->StatusOutlet, "error", "descriptor names must be symbols");
            return;
        }
        object->State->Config.Descriptors.emplace_back(atom_getsym(argv + index)->s_name);
    }
    max_train_apply_configuration(object);
}

static void max_train_set_training_sample_rates(MaxOpenScofoTrain *object, t_symbol *, long argc, t_atom *argv) {
    if (!max_train_can_configure(object))
        return;
    if (argc == 0) {
        max_train_output(object->StatusOutlet, "error", "training_sample_rates requires at least one rate");
        return;
    }
    object->State->Config.TrainingSampleRates.clear();
    for (long index = 0; index < argc; ++index) {
        if (!max_train_is_number(argv[index])) {
            max_train_output(object->StatusOutlet, "error", "training sample rates must be numbers");
            return;
        }
        object->State->Config.TrainingSampleRates.push_back(static_cast<int>(atom_getfloat(argv + index)));
    }
    max_train_apply_configuration(object);
}

#define OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(Name, Field, Type)                                                           \
    static void Name(MaxOpenScofoTrain *object, double value) {                                                        \
        if (!max_train_can_configure(object))                                                                          \
            return;                                                                                                    \
        object->State->Config.Field = static_cast<Type>(value);                                                        \
        max_train_apply_configuration(object);                                                                         \
    }

OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_sample_rate, SampleRate, int)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_fft_size, FFTSize, int)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_hop_size, HopSize, int)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_min_db, MinDB, double)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_max_pitch_confidence, MaxPitchConfidence, double)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_augmentation_factor, AugmentationFactor, double)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_test_fraction, TestFraction, double)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_iterations, Iterations, int)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_learning_rate, LearningRate, double)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_depth, Depth, int)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_early_stopping, EarlyStoppingRounds, int)
OPENSCOFO_MAX_TRAIN_NUMBER_SETTER(max_train_set_seed, RandomSeed, int)

#undef OPENSCOFO_MAX_TRAIN_NUMBER_SETTER

static void max_train_analyze(MaxOpenScofoTrain *object) {
    max_train_start_job(object, MaxTrainJob::Analyze);
}

static void max_train_train(MaxOpenScofoTrain *object) {
    max_train_start_job(object, MaxTrainJob::Train);
}

static void max_train_export(MaxOpenScofoTrain *object, t_symbol *path) {
    if (path == nullptr || path->s_name == nullptr || path->s_name[0] == '\0') {
        max_train_output(object->StatusOutlet, "error", "export requires one model path");
        return;
    }
    max_train_start_job(object, MaxTrainJob::Export, path->s_name);
}

static void max_train_cut(MaxOpenScofoTrain *object, t_symbol *, long argc, t_atom *argv) {
    if (!max_train_can_configure(object))
        return;
    if (argc != 2 || atom_gettype(argv) != A_SYM || atom_gettype(argv + 1) != A_SYM) {
        max_train_output(object->StatusOutlet, "error", "cut requires a buffer and a label");
        return;
    }

    t_symbol *bufferName = atom_getsym(argv);
    t_buffer_ref *bufferReference = buffer_ref_new(reinterpret_cast<t_object *>(object), bufferName);
    if (bufferReference == nullptr) {
        max_train_output(object->StatusOutlet, "error", "failed to create buffer reference");
        return;
    }
    t_buffer_obj *buffer = buffer_ref_getobject(bufferReference);
    if (buffer == nullptr) {
        object_error(reinterpret_cast<t_object *>(object), "buffer %s was not found", bufferName->s_name);
        object_free(bufferReference);
        return;
    }

    const t_atom_long frameCount = buffer_getframecount(buffer);
    const t_atom_long channelCount = buffer_getchannelcount(buffer);
    const double bufferSampleRate = buffer_getsamplerate(buffer);
    if (frameCount <= 0 || channelCount <= 0 || bufferSampleRate <= 0.0 ||
        bufferSampleRate > static_cast<double>(std::numeric_limits<int>::max())) {
        max_train_output(object->StatusOutlet, "error", "buffer has invalid dimensions or sample rate");
        object_free(bufferReference);
        return;
    }
    float *bufferSamples = buffer_locksamples(buffer);
    if (bufferSamples == nullptr) {
        max_train_output(object->StatusOutlet, "error", "failed to lock buffer samples");
        object_free(bufferReference);
        return;
    }

    std::vector<float> samples(static_cast<std::size_t>(frameCount));
    const float channelScale = 1.0F / static_cast<float>(channelCount);
    for (t_atom_long frame = 0; frame < frameCount; ++frame) {
        float sample = 0.0F;
        for (t_atom_long channel = 0; channel < channelCount; ++channel) {
            const std::size_t sampleIndex = static_cast<std::size_t>(frame) * static_cast<std::size_t>(channelCount) +
                                            static_cast<std::size_t>(channel);
            sample += bufferSamples[sampleIndex];
        }
        samples[static_cast<std::size_t>(frame)] = sample * channelScale;
    }
    buffer_unlocksamples(buffer);
    object_free(bufferReference);

    max_train_start_job(object, MaxTrainJob::Cut, atom_getsym(argv + 1)->s_name, std::move(samples),
                        static_cast<int>(bufferSampleRate));
}

static void max_train_clear_data(MaxOpenScofoTrain *object) {
    max_train_start_job(object, MaxTrainJob::ClearData);
}

static void max_train_cancel(MaxOpenScofoTrain *object) {
    if (!object->State->Running.load()) {
        max_train_output(object->StatusOutlet, "status", "idle");
        return;
    }
    object->State->Trainer.Cancel();
    max_train_output(object->StatusOutlet, "status", "cancelling");
}

static void max_train_assist(MaxOpenScofoTrain *, void *, long message, long argument, char *text) {
    if (message == ASSIST_INLET) {
        snprintf(text, 256, "Training configuration and commands");
        return;
    }
    if (argument == 0)
        snprintf(text, 256, "Progress");
    else if (argument == 1)
        snprintf(text, 256, "Status and errors");
    else
        snprintf(text, 256, "Evaluation and segment counts");
}

static void *max_train_new(t_symbol *, long, t_atom *) {
    auto *object = static_cast<MaxOpenScofoTrain *>(object_alloc(max_train_class));
    if (object == nullptr) {
        return nullptr;
    }
    object->State = nullptr;
    object->Clock = nullptr;
    object->ProgressOutlet = nullptr;
    object->StatusOutlet = nullptr;
    object->EvaluationOutlet = nullptr;

    char patchPath[MAX_PATH_CHARS] = {};
    if (path_toabsolutesystempath(path_getdefault(), "", patchPath) != 0) {
        std::snprintf(patchPath, sizeof(patchPath), ".");
    }
    try {
        object->State = new MaxTrainState(patchPath);
    } catch (const std::exception &error) {
        object_error(reinterpret_cast<t_object *>(object), "%s", error.what());
        object_free(object);
        return nullptr;
    }

    object->EvaluationOutlet = outlet_new(object, nullptr);
    object->StatusOutlet = outlet_new(object, nullptr);
    object->ProgressOutlet = outlet_new(object, nullptr);
    object->Clock = clock_new(object, reinterpret_cast<method>(max_train_tick));
    return object;
}

static void max_train_free(MaxOpenScofoTrain *object) {
    if (object->Clock != nullptr) {
        clock_unset(object->Clock);
    }
    delete object->State;
    object->State = nullptr;
    if (object->Clock != nullptr) {
        object_free(object->Clock);
        object->Clock = nullptr;
    }
}

int C74_EXPORT main() {
    t_class *external = class_new("openscofo.train", reinterpret_cast<method>(max_train_new),
                                  reinterpret_cast<method>(max_train_free), sizeof(MaxOpenScofoTrain), 0L, A_GIMME, 0);

    class_addmethod(external, reinterpret_cast<method>(max_train_set_folder), "folder", A_GIMME, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_ir_folders), "ir_folder", A_GIMME, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_descriptors), "descriptors", A_GIMME, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_sample_rate), "sample_rate", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_training_sample_rates), "training_sample_rates",
                    A_GIMME, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_fft_size), "fft_size", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_hop_size), "hop_size", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_min_db), "min_db", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_max_pitch_confidence), "max_pitch_confidence",
                    A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_augmentation_factor), "augmentation_factor",
                    A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_test_fraction), "test_fraction", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_iterations), "iterations", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_learning_rate), "learning_rate", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_depth), "depth", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_early_stopping), "early_stopping_rounds", A_FLOAT,
                    0);
    class_addmethod(external, reinterpret_cast<method>(max_train_set_seed), "seed", A_FLOAT, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_analyze), "analyze", A_NOTHING, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_train), "train", A_NOTHING, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_export), "export", A_SYM, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_cut), "cut", A_GIMME, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_clear_data), "clear_data", A_NOTHING, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_cancel), "cancel", A_NOTHING, 0);
    class_addmethod(external, reinterpret_cast<method>(max_train_assist), "assist", A_CANT, 0);

    class_register(CLASS_BOX, external);
    max_train_class = external;
    object_post(nullptr, "[openscofo.train] native trainer %s", OPENSCOFO_VERSION);
    return 0;
}
