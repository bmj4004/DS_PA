import re
import matplotlib.pyplot as plt
import numpy as np

def parse_data(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    
    memory_usages = re.findall(r'최대: (\d+) KB', content)
    return list(map(int, memory_usages))

# 데이터 파싱
file1_memory = parse_data('Only_3D_Array_mem.txt')
file2_memory = parse_data('2D_Array_Min_Heap_mem.txt')

# 이상치 제거 (Q1 - 1.5 * IQR, Q3 + 1.5 * IQR 범위 밖의 값)
def remove_outliers(data):
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - (1.5 * iqr)
    upper_bound = q3 + (1.5 * iqr)
    return [x for x in data if lower_bound <= x <= upper_bound]

file1_memory = remove_outliers(file1_memory)
file2_memory = remove_outliers(file2_memory)

# 박스플롯 그리기
fig, ax = plt.subplots(figsize=(8, 6))
ax.boxplot([file1_memory, file2_memory], labels=['Only_3D_Array', '2D_Array_Min_Heap'])
ax.set_ylabel('Memory Usage (KB)')
ax.set_title('Comparison of Maximum Memory Usage')

plt.tight_layout()

# PNG로 저장
plt.savefig('memory_usage_comparison.png')

plt.show()