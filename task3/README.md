# Отчет по анализу данных секвенирования E.coli с использованием Flyte

## 1. Ссылка на загруженные прочтения из NCBI SRA
Используется рид-сет: SRR33692911 (пример WSG e.coil)
Ссылка на NCBI SRA: https://www.ncbi.nlm.nih.gov/sra/SRR33692911


## 2. Скрипт на bash с реализованным алгоритмом
**Файл:** `scripts/analysis.sh`

```bash
#!/bin/bash

# Параметры
REF="data/ecoli_ref.fna"
FASTQ="data/SRR33692911.fastq"

# Контроль качества
fastqc $FASTQ -o fastqc_results

# Выравнивание
bwa mem -t 4 $REF $FASTQ > aligned.sam

# Статистика
samtools flagstat aligned.sam > flagstat.txt

# Анализ
python3 scripts/parse_flagstat.py flagstat.txt
```

## 3. Результат команды samtools flagstat
**Файл:** `results/flagstat.txt`
```
1165545 + 0 in total (QC-passed reads + QC-failed reads)
910504 + 0 mapped (78.12% : N/A)
903417 + 0 primary mapped (77.98% : N/A)
... [полный вывод]
```

## 4. Скрипт разбора результатов
**Файл:** `scripts/parse_flagstat.py`

## 6. Инструкция по установке Flyte
**Установка:**
```bash
# 1. Установите Docker
sudo apt install docker.io
sudo systemctl enable --now docker

# 2. Установите flytectl
curl -sL https://raw.githubusercontent.com/flyteorg/flytectl/master/install.sh | bash
sudo mv ~/.flyte/bin/flytectl /usr/local/bin/

# 3. Запустите кластер
flytectl demo start
```

## 7. Тестовый пайплайн ("Hello World")
**Файл:** `flyte/hello_world.py`
```python
from flytekit import task, workflow

@task
def say_hello(name: str) -> str:
    return f"Hello, {name}!"

@workflow
def hello_wf(name: str = "World") -> str:
    return say_hello(name=name)
```

## 8. Результаты работы пайплайна
**Вывод:**
```
Hello, Bioinformatics!
```

## 9. Инструменты визуализации
С визуализацией не получилось у меня ничего

## 10. Пайплайн оценки качества
**Файл:** `flyte/ecoli_analysis.py`  

## 11. Результаты работы пайплайна
```
Mapping percentage: 78.12%
Status: PASS
```
