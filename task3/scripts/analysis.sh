#!/bin/bash

# Параметры
REF="data/ecoli_ref.fna"
FASTQ="data/SRR33692911.fastq"
THRESHOLD=75

# 1. Проверка наличия инструментов
command -v fastqc >/dev/null 2>&1 || { echo >&2 "FastQC не установлен. Установите сначала."; exit 1; }
command -v bwa >/dev/null 2>&1 || { echo >&2 "BWA не установлен. Установите сначала."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "Samtools не установлен. Установите сначала."; exit 1; }

# 2. Контроль качества
echo "Запуск FastQC..."
mkdir -p fastqc_results
fastqc "$FASTQ" -o fastqc_results || { echo "Ошибка FastQC"; exit 1; }

# 3. Индексация референса (если ещё не сделано)
if [ ! -f "${REF}.bwt" ]; then
    echo "Индексация референсного генома..."
    bwa index "$REF" || { echo "Ошибка индексации"; exit 1; }
fi

# 4. Выравнивание
echo "Выравнивание reads..."
bwa mem -t 4 "$REF" "$FASTQ" > aligned.sam || { echo "Ошибка выравнивания"; exit 1; }

# 5. Конвертация и сортировка
echo "Конвертация SAM в BAM..."
samtools view -@ 4 -Sb aligned.sam > aligned.bam || { echo "Ошибка конвертации"; exit 1; }

echo "Сортировка BAM..."
samtools sort -@ 4 -o aligned_sorted.bam aligned.bam || { echo "Ошибка сортировки"; exit 1; }

# 6. Статистика
echo "Расчет статистики..."
samtools flagstat aligned_sorted.bam > flagstat.txt || { echo "Ошибка flagstat"; exit 1; }

# 7. Анализ результатов
echo "Анализ результатов..."
python3 scripts/parse_flagstat.py flagstat.txt || { echo "Ошибка анализа"; exit 1; }

echo "Анализ завершен успешно!"
