import sys
import re

def parse_flagstat(file_path):
    with open(file_path) as f:
        data = f.read()
    
    # Ищем строку с процентами картирования
    mapped_line = next(line for line in data.split('\n') if "mapped (" in line)
    
    # Извлекаем процент
    percent = float(mapped_line.split("(")[1].split("%")[0])
    total = int(mapped_line.split()[0])
    mapped = int(mapped_line.split()[0])  # Первое число в строке
    
    print(f"Результаты анализа:")
    print(f"Всего reads: {total}")
    print(f"Картировано: {mapped} ({percent}%)")
    return percent

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Использование: python parse_flagstat.py <flagstat.txt>")
        sys.exit(1)
    
    try:
        percent = parse_flagstat(sys.argv[1])
        sys.exit(0 if percent >= 75 else 1)
    except Exception as e:
        print(f"Ошибка: {str(e)}")
        sys.exit(1)

