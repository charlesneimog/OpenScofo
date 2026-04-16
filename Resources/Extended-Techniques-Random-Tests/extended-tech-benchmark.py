#!/usr/bin/env python3
"""
Extended Techniques Benchmark for OpenScofo
Validação de implementações com persistência e comparação automática
"""

import os
import json
import hashlib
from datetime import datetime
from pathlib import Path
from typing import List, Tuple, Dict
import numpy as np
from scipy import stats
import librosa
import OpenScofo

os.chdir(os.path.dirname(__file__))

# ============================================================================
# CONFIGURAÇÃO GLOBAL
# ============================================================================
SR = 48000
FFT = 2048
HOP = 512

# Caminhos dos arquivos de teste
TEST_FILES = [
    {"audio": "./score-1.wav", "score": "./score-1.txt"},
    {"audio": "./score-2.wav", "score": "./score-2.txt"},
    {"audio": "./score-3.wav", "score": "./score-3.txt"},
]

RESULTS_PATH = "follower_validation.json"
EPSILON = 1e-6  # Tolerância para considerar scores iguais


# ============================================================================
# CLASSE VALIDADORA
# ============================================================================
class ScoreFollowerValidator:
    def __init__(self, results_path: str = RESULTS_PATH):
        self.results_path = Path(results_path)
        self.current_results = None
        self.previous_results = None

    def compute_metrics(self, events: List[Tuple]) -> Dict:
        """Compute comprehensive metrics beyond simple mean/std"""
        if not events:
            return {"error": "no_events", "composite_score": float("inf")}

        errors = np.array([e[3] for e in events])  # em segundos
        abs_errors = np.abs(errors)

        # Métricas básicas
        metrics = {
            "n_events": len(events),
            "mean_error_ms": float(np.mean(errors) * 1000),
            "std_error_ms": float(np.std(errors) * 1000),
            "mae_ms": float(np.mean(abs_errors) * 1000),
            "median_error_ms": float(np.median(errors) * 1000),
            "median_absolute_error_ms": float(np.median(abs_errors) * 1000),
            "rmse_ms": float(np.sqrt(np.mean(errors**2)) * 1000),
            "percentile_95_abs_ms": float(np.percentile(abs_errors, 95) * 1000),
            "percentile_99_abs_ms": float(np.percentile(abs_errors, 99) * 1000),
        }

        # Assimetria do erro
        negative_errors = errors[errors < 0]
        positive_errors = errors[errors > 0]

        metrics["bias_ms"] = float(np.mean(errors) * 1000)
        metrics["negative_bias_ms"] = (
            float(np.mean(negative_errors) * 1000) if len(negative_errors) > 0 else 0
        )
        metrics["positive_bias_ms"] = (
            float(np.mean(positive_errors) * 1000) if len(positive_errors) > 0 else 0
        )
        metrics["ratio_negative"] = float(len(negative_errors) / len(errors))

        # Detecção de outliers catastróficos
        if len(errors) > 3:
            q75, q25 = np.percentile(abs_errors, 75), np.percentile(abs_errors, 25)
            iqr = q75 - q25
            if iqr > 0:
                outliers = abs_errors > (q75 + 1.5 * iqr)
                metrics["outlier_ratio"] = float(np.sum(outliers) / len(errors))
            else:
                metrics["outlier_ratio"] = 0.0
        else:
            metrics["outlier_ratio"] = 0.0

        # Score composto (menor = melhor)
        normalized_mae = min(metrics["mae_ms"] / 100.0, 2.0)
        normalized_outliers = min(metrics["outlier_ratio"] * 10, 2.0)
        normalized_bias = min(abs(metrics["bias_ms"]) / 50.0, 1.0)

        metrics["composite_score"] = (
            0.5 * normalized_mae + 0.3 * normalized_outliers + 0.2 * normalized_bias
        )

        return metrics

    def hash_implementation(self) -> str:
        """Create hash of the compiled OpenScofo module"""
        try:
            import OpenScofo

            module_path = Path(OpenScofo.__file__)
            if module_path.exists():
                with open(module_path, "rb") as f:
                    content = f.read()
                    return hashlib.sha256(content).hexdigest()[:8]
        except:
            pass

        # Fallback: timestamp da última execução
        return hashlib.sha256(str(datetime.now()).encode()).hexdigest()[:8]

    def save_results(
        self,
        events: List[Tuple],
        audio_file: str,
        score_file: str,
        implementation_name: str = None,
    ):
        """Save results with comparison logic"""
        metrics = self.compute_metrics(events)

        if implementation_name is None:
            impl_hash = self.hash_implementation()
            implementation_name = f"impl_{impl_hash}"

        result_entry = {
            "timestamp": datetime.now().isoformat(),
            "implementation": implementation_name,
            "audio_file": audio_file,
            "score_file": score_file,
            "metrics": metrics,
            "n_events": len(events),
        }

        # Carregar resultados anteriores
        if self.results_path.exists():
            with open(self.results_path, "r") as f:
                previous_data = json.load(f)
        else:
            previous_data = {"history": [], "best_implementations": {}}

        # Adicionar resultado atual
        previous_data["history"].append(result_entry)

        # Rastrear melhor implementação para este arquivo específico
        if audio_file not in previous_data["best_implementations"]:
            previous_data["best_implementations"][audio_file] = implementation_name
            is_better = True
            is_equal = False
            comparison_msg = f"📊 Primeira execução para {audio_file} (baseline)"
        else:
            best_impl = previous_data["best_implementations"][audio_file]
            best_metrics = None

            for entry in previous_data["history"]:
                if (
                    entry["implementation"] == best_impl
                    and entry["audio_file"] == audio_file
                ):
                    best_metrics = entry["metrics"]
                    break

            if best_metrics:
                score_diff = (
                    metrics["composite_score"] - best_metrics["composite_score"]
                )

                if abs(score_diff) < EPSILON:
                    # Praticamente igual
                    is_better = False
                    is_equal = True
                    comparison_msg = f"≈ IGUAL (score: {metrics['composite_score']:.3f} vs {best_metrics['composite_score']:.3f})"
                elif metrics["composite_score"] < best_metrics["composite_score"]:
                    is_better = True
                    is_equal = False
                    previous_data["best_implementations"][
                        audio_file
                    ] = implementation_name
                    comparison_msg = f"✓ MELHOR (score: {metrics['composite_score']:.3f} vs {best_metrics['composite_score']:.3f} [-{score_diff:.3f}])"
                else:
                    is_better = False
                    is_equal = False
                    comparison_msg = f"✗ PIOR (score: {metrics['composite_score']:.3f} vs {best_metrics['composite_score']:.3f} [+{abs(score_diff):.3f}])"
            else:
                is_better = True
                is_equal = False
                comparison_msg = f"📊 Primeira referência para {audio_file}"

            print(f"\n--- COMPARAÇÃO para {Path(audio_file).name} ---")
            print(comparison_msg)
            if best_metrics:
                print(
                    f"  MAE: {metrics['mae_ms']:.2f}ms vs {best_metrics['mae_ms']:.2f}ms"
                )
                print(
                    f"  Outliers: {metrics['outlier_ratio']:.2%} vs {best_metrics['outlier_ratio']:.2%}"
                )
                print(
                    f"  Bias: {metrics['bias_ms']:+.2f}ms vs {best_metrics['bias_ms']:+.2f}ms"
                )

        # Salvar
        with open(self.results_path, "w") as f:
            json.dump(previous_data, f, indent=2)

        print(f"\nResultados salvos em {self.results_path}")
        return is_better


# ============================================================================
# FUNÇÃO PRINCIPAL DE PROCESSAMENTO
# ============================================================================
def process_audio_file(audio_path: str, score_path: str) -> List[Tuple]:
    """
    Processa um arquivo de áudio com o OpenScofo e retorna eventos detectados
    Returns: List of (score_pos, detected_time, expected_time, error)
    """
    print(f"\n--- Processando {audio_path} ---")

    # Load audio
    x, sr = librosa.load(audio_path, sr=SR)

    # Init Scofo
    scofo = OpenScofo.OpenScofo(SR, FFT, HOP)
    scofo.parse_score(Path(score_path))

    states = scofo.get_states()

    # Map score_pos -> expected onset time
    score_pos_to_expected = {}
    for s in states:
        if hasattr(s, "score_pos") and hasattr(s, "onset_expected"):
            score_pos_to_expected[s.score_pos] = s.onset_expected

    # Process
    n = len(x)
    prev_pos = None
    events = []

    for frame_idx, start in enumerate(range(0, n, HOP)):
        end = start + HOP

        if end <= n:
            frame = x[start:end]
        else:
            frame = np.zeros(FFT, dtype=x.dtype)
            valid = n - start
            if valid > 0:
                frame[:valid] = x[start:n]

        scofo.process_block(frame)
        pos = scofo.get_current_score_position()

        if pos < 0:
            continue

        if pos != prev_pos:
            t_det = start / SR
            t_exp = score_pos_to_expected.get(pos, None)

            if t_exp is not None:
                err = t_det - t_exp
                events.append((pos, t_det, t_exp, err))

                print(
                    f"pos={pos:03d}  det={t_det:8.3f}s  exp={t_exp:8.3f}s  err={err*1000:+8.2f} ms"
                )

            prev_pos = pos

    print(f"Total eventos detectados: {len(events)}")
    return events


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 60)
    print("EXTENDED TECHNIQUES BENCHMARK - OpenScofo")
    print("=" * 60)

    # Mostrar versão do módulo compilado
    try:
        if hasattr(OpenScofo, "__file__"):
            module_time = Path(OpenScofo.__file__).stat().st_mtime
            print(f"Módulo OpenScofo: {Path(OpenScofo.__file__).name}")
            print(
                f"Última modificação: {datetime.fromtimestamp(module_time).strftime('%Y-%m-%d %H:%M:%S')}"
            )
    except:
        pass
    print()

    validator = ScoreFollowerValidator(RESULTS_PATH)

    # Processar cada arquivo de teste
    all_results = []

    for test_file in TEST_FILES:
        audio_path = test_file["audio"]
        score_path = test_file["score"]

        # Verificar se arquivos existem
        if not Path(audio_path).exists():
            print(f"⚠️ Arquivo não encontrado: {audio_path}")
            continue

        if not Path(score_path).exists():
            print(f"⚠️ Arquivo não encontrado: {score_path}")
            continue

        # Processar
        events = process_audio_file(audio_path, score_path)

        if events:
            is_better = validator.save_results(events, audio_path, score_path)
            all_results.append(
                {"file": audio_path, "is_better": is_better, "n_events": len(events)}
            )
        else:
            print(f"❌ Nenhum evento detectado em {audio_path}")

    # Sumário final
    print("\n" + "=" * 60)
    print("RESUMO DA VALIDAÇÃO")
    print("=" * 60)

    if not all_results:
        print("❌ Nenhum arquivo processado com sucesso")
        return

    better_count = sum(1 for r in all_results if r["is_better"])
    total = len(all_results)

    print(f"Arquivos testados: {total}")
    print(f"Melhorias detectadas: {better_count}/{total}")
    print(f"Piores ou iguais: {total - better_count}/{total}")

    if better_count == total:
        print("\n✅ IMPLEMENTAÇÃO APROVADA - Melhor em TODOS os testes")
    elif better_count > total / 2:
        print("\n✅ IMPLEMENTAÇÃO APROVADA - Melhor na maioria dos testes")
    elif better_count == 0:
        print("\n≈ IMPLEMENTAÇÃO EQUIVALENTE - Sem diferença significativa")
        print("   As mudanças não afetaram a performance (positiva ou negativamente)")
    else:
        print("\n⚠️ RESULTADO INCONCLUSIVO - Mais testes necessários")


if __name__ == "__main__":
    main()
