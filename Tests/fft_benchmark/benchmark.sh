#!/bin/bash

# Configurações do Teste
RUNS=20
BIN_FFTW3="./build/bench_fftw3"
BIN_PFFFT="./build/bench_pffft"

echo "==================================================="
echo " Iniciando Bateria de Testes ($RUNS iterações)"
echo " Cada iteração executa 1.000.000 de FFTs internas."
echo "==================================================="

# Função para executar os binários e extrair a média
run_and_parse() {
    local binary=$1
    local name=$2
    local total=0
    
    echo -n "Avaliando $name "
    
    for i in $(seq 1 $RUNS); do
        # Executa o binário, pega a linha da média e extrai o 6º bloco de texto (o número)
        val=$($binary | grep "Average Time" | awk '{print $6}')
        
        # Soma o valor ao total usando a calculadora de precisão do bash (bc)
        total=$(echo "$total + $val" | bc -l)
        
        # Imprime um ponto para mostrar progresso
        echo -n "."
    done
    
    # Calcula a média real dividindo o total pelo número de execuções
    local avg=$(echo "$total / $RUNS" | bc -l)
    
    # Formata a saída para 5 casas decimais
    printf "\n-> Média Estável do %s: %.5f us\n\n" "$name" "$avg"
}

# Executa as funções
run_and_parse $BIN_FFTW3 "FFTW3"
run_and_parse $BIN_PFFFT "PFFFT"

echo "==================================================="
echo " Bateria concluída."
echo "==================================================="
