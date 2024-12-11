import re
import matplotlib.pyplot as plt
import numpy as np

def parse_data(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    
    grid_times = re.findall(r'그리드 구축: (\d+\.\d+) ms', content)
    png_times = re.findall(r'PNG 생성: (\d+\.\d+) ms', content)
    total_times = re.findall(r'전체 처리: (\d+\.\d+) ms', content)
    
    return list(map(float, grid_times)), list(map(float, png_times)), list(map(float, total_times))

# 데이터 파싱
only_3d_array_grid, only_3d_array_png, only_3d_array_total = parse_data('Only_3D_Array.txt')
min_heap_grid, min_heap_png, min_heap_total = parse_data('2D_Array_Min_Heap.txt')

# 이상치 제거 (Q1 - 1.5 * IQR, Q3 + 1.5 * IQR 범위 밖의 값)
def remove_outliers(data):
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - (1.5 * iqr)
    upper_bound = q3 + (1.5 * iqr)
    return [x for x in data if lower_bound <= x <= upper_bound]

only_3d_array_grid = remove_outliers(only_3d_array_grid)
only_3d_array_png = remove_outliers(only_3d_array_png)  
only_3d_array_total = remove_outliers(only_3d_array_total)
min_heap_grid = remove_outliers(min_heap_grid)
min_heap_png = remove_outliers(min_heap_png)
min_heap_total = remove_outliers(min_heap_total)

# 박스플롯 그리기 및 PNG로 저장
fig, ax = plt.subplots(figsize=(6, 6))
bp = ax.boxplot([only_3d_array_grid, min_heap_grid], widths=0.6)
ax.set_xticklabels(['Only_3D_Array', '2D_Array_Min_Heap'])
ax.set_ylabel('Time (ms)')
ax.set_title('Grid Construction')
for median in bp['medians']:
    median.set(color='red', linewidth=2)
plt.tight_layout()
plt.savefig('grid_construction.png')
plt.close()

fig, ax = plt.subplots(figsize=(6, 6))
bp = ax.boxplot([only_3d_array_png, min_heap_png], widths=0.6)
ax.set_xticklabels(['Only_3D_Array', '2D_Array_Min_Heap'])
ax.set_ylabel('Time (ms)')
ax.set_title('2D Projection')
for median in bp['medians']:
    median.set(color='red', linewidth=2)
plt.tight_layout()
plt.savefig('png_generation.png')
plt.close()

fig, ax = plt.subplots(figsize=(6, 6))
bp = ax.boxplot([only_3d_array_total, min_heap_total], widths=0.6)
ax.set_xticklabels(['Only_3D_Array', '2D_Array_Min_Heap'])  
ax.set_ylabel('Time (ms)')
ax.set_title('Total Processing')
for median in bp['medians']:
    median.set(color='red', linewidth=2)
plt.tight_layout() 
plt.savefig('total_processing.png')
plt.close()