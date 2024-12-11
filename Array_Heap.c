#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <png.h>
#include <libgen.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <malloc.h>

// 설정 가능한 파라미터
#define HORIZONTAL_FOV 70.0f
#define VERTICAL_FOV 43.0f
#define ANGULAR_RESOLUTION 0.5f
#define CENTER_ELEVATION 90.0f
#define CAMERA_HEIGHT 1.0f
#define NUM_ITERATIONS 100

// 전역 변수로 메모리 추적
static long current_memory = 0;  // 현재 할당된 메모리 (bytes)
static long peak_memory = 0;     // 최대 메모리 사용량 (bytes)

// 매크로 정의 - 노드와 그 자식들의 평균 거리 계산
#define CALC_AVG_DISTANCE(node) \
        ((node)->distance + \
        ((node)->left ? (node)->left->distance : (node)->distance) + \
        ((node)->right ? (node)->right->distance : (node)->distance)) / 3.0f

// 성능 측정을 위한 구조체
typedef struct {
    double time_taken;    // 소요 시간
    long memory_used;     // 사용된 메모리
} PerfMetrics;

typedef struct {
    PerfMetrics grid_construction;
    PerfMetrics png_generation;
    PerfMetrics total;
    long baseline_memory;
} ProcessingMetrics;

// 색상 배열 정의 (가까운 순서대로, 15m까지)
const float colors[][3] = {
    {1.00f, 0.00f, 0.00f},  // 0-1m  (빨강)
    {1.00f, 0.25f, 0.00f},  // 1-2m  (주황빨강)
    {1.00f, 0.50f, 0.00f},  // 2-3m  (주황)
    {1.00f, 0.75f, 0.00f},  // 3-4m  (주황노랑)
    {1.00f, 1.00f, 0.00f},  // 4-5m  (노랑)
    {0.75f, 1.00f, 0.00f},  // 5-6m  (연두)
    {0.50f, 1.00f, 0.00f},  // 6-7m  (연두초록)
    {0.00f, 1.00f, 0.00f},  // 7-8m  (초록)
    {0.00f, 1.00f, 0.50f},  // 8-9m  (청록)
    {0.00f, 1.00f, 0.75f},  // 9-10m (밝은청록)
    {0.00f, 1.00f, 1.00f},  // 10-11m (하늘)
    {0.00f, 0.75f, 1.00f},  // 11-12m (밝은파랑)
    {0.00f, 0.50f, 1.00f},  // 12-13m (파랑)
    {0.00f, 0.25f, 1.00f},  // 13-14m (진한파랑)
    {0.00f, 0.00f, 1.00f},  // 14-15m (더진한파랑)
    {0.00f, 0.00f, 0.75f}   // 15m 이상 (남색)
};

// 힙 노드 구조체
typedef struct HeapNode {
    float x, y, z;            // 직교 좌표
    float distance;           // 원점 거리
    float azimuth;           // 방위각 (0~360도)
    float elevation;         // 고도각 (0~180도)
    struct HeapNode* left;    // 왼쪽 자식
    struct HeapNode* right;   // 오른쪽 자식
} HeapNode;

// Min Heap 구조체
typedef struct {
    HeapNode* root;    // 루트 노드
    int size;          // 힙의 크기
} MinHeap;

// 구면 그리드 셀 구조체
typedef struct {
    MinHeap heap;      // min heap
    int is_valid;      // FOV 내부 여부
} SphericalCell;

// 구면 그리드 구조체
typedef struct {
    SphericalCell** cells;  // 2D 배열
    int azimuth_cells;    
    int elevation_cells;  
} SphericalGrid;

// malloc 래퍼 함수
void* tracked_malloc(size_t size) {
    void* ptr = malloc(size);
    if (ptr) {
        current_memory += size;
        if (current_memory > peak_memory) {
            peak_memory = current_memory;
        }
    }
    return ptr;
}

// free 래퍼 함수
void tracked_free(void* ptr, size_t size) {
    if (ptr) {
        current_memory -= size;
        free(ptr);
    }
}

// 메모리 사용량 조회 함수
long get_memory_usage() {
    return current_memory / 1024;  // bytes를 KB로 변환하여 반환
}

// 메모리 추적 초기화
void reset_memory_tracking() {
    current_memory = 0;
    peak_memory = 0;
}

// 현재 시간을 마이크로초 정밀도로 반환
double get_current_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

// 표준편차 계산
double calculate_std_dev(double values[], int n, double mean) {
    double sum_sq_diff = 0;
    for (int i = 0; i < n; i++) {
        double diff = values[i] - mean;
        sum_sq_diff += diff * diff;
    }
    return sqrt(sum_sq_diff / n);
}

// 통계 출력 (파일과 터미널 모두 지원)
void print_statistics(const char* phase, PerfMetrics metrics[], int n, FILE* file, int print_to_terminal) {
    double times[NUM_ITERATIONS];
    double memories[NUM_ITERATIONS];
    double time_sum = 0, memory_sum = 0;
    double time_max = metrics[0].time_taken * 1000.0;
    double memory_max = metrics[0].memory_used;
    
    for (int i = 0; i < n; i++) {
        times[i] = metrics[i].time_taken * 1000.0;
        memories[i] = metrics[i].memory_used;
        time_sum += times[i];
        memory_sum += memories[i];
        
        time_max = fmax(time_max, times[i]);
        memory_max = fmax(memory_max, memories[i]);
    }
    
    double time_mean = time_sum / n;
    double memory_mean = memory_sum / n;
    double time_std = calculate_std_dev(times, n, time_mean);
    double memory_std = calculate_std_dev(memories, n, memory_mean);
    
    if (print_to_terminal) {
        printf("\n=== %s 통계 (n=%d) ===\n", phase, n);
        printf("처리 시간:\n");
        printf("  평균: %.2f ms ± %.2f ms\n", time_mean, time_std);
        printf("  최대: %.2f ms\n", time_max);
        printf("메모리 사용량:\n");
        printf("  평균: %ld KB ± %.2f KB\n", (long)memory_mean, memory_std);
        printf("  최대: %ld KB\n", (long)memory_max);
    }
    
    if (file) {
        fprintf(file, "\n=== %s 통계 (n=%d) ===\n", phase, n);
        fprintf(file, "처리 시간:\n");
        fprintf(file, "  평균: %.2f ms ± %.2f ms\n", time_mean, time_std);
        fprintf(file, "  최대: %.2f ms\n", time_max);
        fprintf(file, "메모리 사용량:\n");
        fprintf(file, "  평균: %ld KB ± %.2f KB\n", (long)memory_mean, memory_std);
        fprintf(file, "  최대: %ld KB\n", (long)memory_max);
    }
}

// float 값을 리틀 엔디안으로 읽는 함수
float read_float_le(FILE* file) {
    unsigned char buffer[4];
    if (fread(buffer, 1, 4, file) != 4) {
        return 0.0f;
    }
    unsigned int val = ((unsigned int)buffer[3] << 24) |
                      ((unsigned int)buffer[2] << 16) |
                      ((unsigned int)buffer[1] << 8) |
                      buffer[0];
    return *(float*)&val;
}

// 새로운 노드 생성
HeapNode* create_node(float x, float y, float z) {
    HeapNode* node = (HeapNode*)tracked_malloc(sizeof(HeapNode));
    if (!node) return NULL;
    
    node->x = x;
    node->y = y;
    node->z = z;
    node->distance = sqrt(x*x + y*y + z*z);
    node->azimuth = atan2(y, x) * 180.0f / M_PI;
    if(node->azimuth < 0) node->azimuth += 360.0f;
    node->elevation = acos(z / node->distance) * 180.0f / M_PI;
    
    node->left = NULL;
    node->right = NULL;
    
    return node;
}

// 힙에 새로운 노드 삽입을 위한 재귀 함수
HeapNode* insert_recursive(HeapNode* root, HeapNode* new_node, int* level_count, int target_level) {
    if (*level_count > target_level) return root;
    
    if (!root) {
        (*level_count)++;
        return new_node;
    }

    if (root->distance > new_node->distance) {
        // 값 교환
        float temp_x = root->x;
        float temp_y = root->y;
        float temp_z = root->z;
        float temp_dist = root->distance;
        float temp_azimuth = root->azimuth;
        float temp_elevation = root->elevation;

        root->x = new_node->x;
        root->y = new_node->y;
        root->z = new_node->z;
        root->distance = new_node->distance;
        root->azimuth = new_node->azimuth;
        root->elevation = new_node->elevation;

        new_node->x = temp_x;
        new_node->y = temp_y;
        new_node->z = temp_z;
        new_node->distance = temp_dist;
        new_node->azimuth = temp_azimuth;
        new_node->elevation = temp_elevation;
    }

    if (!root->left || (root->right && 
        CALC_AVG_DISTANCE(root->left) > CALC_AVG_DISTANCE(root->right))) {
        root->left = insert_recursive(root->left, new_node, level_count, target_level);
    }
    else {
        root->right = insert_recursive(root->right, new_node, level_count, target_level);
    }

    return root;
}

// 힙에 새로운 노드 삽입
void insert_node(MinHeap* heap, HeapNode* node) {
    int level_count = 0;
    int target_level = (int)log2(heap->size + 1);
    heap->root = insert_recursive(heap->root, node, &level_count, target_level);
    heap->size++;
}

// 구면 그리드 생성
SphericalGrid* create_grid() {
    SphericalGrid* grid = (SphericalGrid*)tracked_malloc(sizeof(SphericalGrid));
    if (!grid) return NULL;
    
    grid->azimuth_cells = (int)(360.0f / ANGULAR_RESOLUTION);
    grid->elevation_cells = (int)(180.0f / ANGULAR_RESOLUTION);
    
    grid->cells = (SphericalCell**)tracked_malloc(grid->azimuth_cells * sizeof(SphericalCell*));
    if (!grid->cells) {
        tracked_free(grid, sizeof(SphericalGrid));
        return NULL;
    }

    for (int i = 0; i < grid->azimuth_cells; i++) {
        grid->cells[i] = (SphericalCell*)tracked_malloc(grid->elevation_cells * sizeof(SphericalCell));
        if (!grid->cells[i]) {
            for (int j = 0; j < i; j++) {
                tracked_free(grid->cells[j], grid->elevation_cells * sizeof(SphericalCell));
            }
            tracked_free(grid->cells, grid->azimuth_cells * sizeof(SphericalCell*));
            tracked_free(grid, sizeof(SphericalGrid));
            return NULL;
        }
        for (int j = 0; j < grid->elevation_cells; j++) {
            grid->cells[i][j].heap.root = NULL;
            grid->cells[i][j].heap.size = 0;
            grid->cells[i][j].is_valid = 0;
        }
    }

    // FOV 영역 설정
    float half_hfov = HORIZONTAL_FOV / 2.0f;
    float half_vfov = VERTICAL_FOV / 2.0f;

    for (int i = 0; i < grid->azimuth_cells; i++) {
        float azimuth = i * ANGULAR_RESOLUTION;
        if (azimuth <= half_hfov || azimuth >= (360.0f - half_hfov)) {
            for (int j = 0; j < grid->elevation_cells; j++) {
                float elevation = j * ANGULAR_RESOLUTION;
                if (fabs(elevation - CENTER_ELEVATION) <= half_vfov) {
                    grid->cells[i][j].is_valid = 1;
                }
            }
        }
    }
    
    return grid;
}

// 노드 추가 함수
void add_point_to_grid(SphericalGrid* grid, float x, float y, float z) {
    float rel_z = z - CAMERA_HEIGHT;
    
    HeapNode* node = create_node(x, y, z);
    if (!node) return;
    
    float distance = sqrt(x*x + y*y + rel_z*rel_z);
    float azimuth = atan2(y, x) * 180.0f / M_PI;
    if(azimuth < 0) azimuth += 360.0f;
    float elevation = acos(rel_z / distance) * 180.0f / M_PI;
    
    if (azimuth < 0 || azimuth >= 360.0f ||
        elevation < 0 || elevation >= 180.0f) {
        tracked_free(node, sizeof(HeapNode));
        return;
    }
    
    int azimuth_idx = (int)(azimuth / ANGULAR_RESOLUTION);
    int elevation_idx = (int)(elevation / ANGULAR_RESOLUTION);
    
    if (azimuth_idx < 0 || azimuth_idx >= grid->azimuth_cells ||
        elevation_idx < 0 || elevation_idx >= grid->elevation_cells) {
        tracked_free(node, sizeof(HeapNode));
        return;
    }

    float half_hfov = HORIZONTAL_FOV / 2.0f;
    float half_vfov = VERTICAL_FOV / 2.0f;

    if ((azimuth <= half_hfov || azimuth >= (360.0f - half_hfov)) &&
        (fabs(elevation - CENTER_ELEVATION) <= half_vfov)) {
        
        if (grid->cells[azimuth_idx][elevation_idx].is_valid) {
            node->distance = distance;
            node->azimuth = azimuth;
            node->elevation = elevation;
            insert_node(&grid->cells[azimuth_idx][elevation_idx].heap, node);
        } else {
            tracked_free(node, sizeof(HeapNode));
        }
    } else {
        tracked_free(node, sizeof(HeapNode));
    }
}

// 재귀적으로 힙 노드들을 해제
void tracked_free_heap_nodes(HeapNode* node) {
    if (node) {
        tracked_free_heap_nodes(node->left);
        tracked_free_heap_nodes(node->right);
        tracked_free(node, sizeof(HeapNode));
    }
}

// 그리드 메모리 해제
void tracked_free_grid(SphericalGrid* grid) {
    if (!grid) return;

    for (int i = 0; i < grid->azimuth_cells; i++) {
        if (grid->cells[i]) {
            for (int j = 0; j < grid->elevation_cells; j++) {
                tracked_free_heap_nodes(grid->cells[i][j].heap.root);
            }
            tracked_free(grid->cells[i], grid->elevation_cells * sizeof(SphericalCell));
        }
    }
    tracked_free(grid->cells, grid->azimuth_cells * sizeof(SphericalCell*));
    tracked_free(grid, sizeof(SphericalGrid));
}

// PCD 파일 처리
SphericalGrid* process_pcd_file(const char* filename) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        printf("파일을 열 수 없습니다: %s\n", filename);
        return NULL;
    }
    
    SphericalGrid* grid = create_grid();
    if (!grid) {
        fclose(file);
        return NULL;
    }
    
    // 헤더 파싱
    char line[1024];
    int num_points = 0;
    while (fgets(line, sizeof(line), file)) {
        if (strncmp(line, "POINTS", 6) == 0) {
            sscanf(line, "POINTS %d", &num_points);
        }
        else if (strcmp(line, "DATA binary\n") == 0) {
            break;
        }
    }

    int valid_points = 0;
    for (int i = 0; i < num_points; i++) {
        float x = read_float_le(file);
        float y = read_float_le(file);
        float z = read_float_le(file);
        unsigned char intensity;
        fread(&intensity, 1, 1, file);
        
        if ((fabs(x) > 0.001f || fabs(y) > 0.001f || fabs(z) > 0.001f) &&
            fabs(x) < 100.0f && fabs(y) < 100.0f && fabs(z) < 100.0f) {
            add_point_to_grid(grid, x, y, z);
            valid_points++;
        }
    }
    
    fclose(file);
    return grid;
}

void save_png(const char* filename, unsigned char* image, int width, int height) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        printf("PNG 파일을 열 수 없습니다: %s\n", filename);
        return;
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fclose(fp);
        return;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        return;
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return;
    }

    png_init_io(png, fp);
    png_set_IHDR(png, info, width, height, 8, PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    png_bytep* rows = (png_bytep*)tracked_malloc(height * sizeof(png_bytep));
    for (int y = 0; y < height; y++) {
        rows[height - 1 - y] = &image[y * width * 3];
    }

    png_write_image(png, rows);
    png_write_end(png, NULL);

    tracked_free(rows, height * sizeof(png_bytep));
    png_destroy_write_struct(&png, &info);
    fclose(fp);
}

// 단일 반복 처리
ProcessingMetrics process_iteration(const char* filename) {
    ProcessingMetrics metrics = {0};
    reset_memory_tracking();  // 메모리 추적 초기화
    
    double start_time = get_current_time();
    long start_memory = current_memory;
    
    // 그리드 구축 단계
    double grid_start = get_current_time();
    SphericalGrid* grid = process_pcd_file(filename);
    double grid_end = get_current_time();
    
    metrics.grid_construction.time_taken = grid_end - grid_start;
    metrics.grid_construction.memory_used = (current_memory - start_memory) / 1024;  // bytes to KB

    if (!grid) {
        printf("PCD 파일 처리 실패\n");
        return metrics;
    }
    
    // PNG 생성 단계
    double png_start = get_current_time();
    long before_png_memory = current_memory;
    
    const int width = 1920;
    const int height = 1080;
    unsigned char* image = (unsigned char*)tracked_malloc(width * height * 3);
    memset(image, 0, width * height * 3);

    // FOV 내 점들만 순회하여 이미지 생성
    for (int i = 0; i < grid->azimuth_cells; i++) {
        for (int j = 0; j < grid->elevation_cells; j++) {
            if (grid->cells[i][j].is_valid && 
                grid->cells[i][j].heap.root) {
                
                HeapNode* node = grid->cells[i][j].heap.root;
                float azimuth = node->azimuth;
                
                float mapped_x;
                if (azimuth >= 324.0f) {
                    mapped_x = (azimuth - 324.0f) / (36.0f) * (width/2);
                } else {
                    mapped_x = (width/2) + (azimuth / 36.0f) * (width/2);
                }
                
                float mapped_y = ((111.0f - node->elevation) / (111.0f - 69.0f)) * height;
                
                int x = width - 1 - (int)mapped_x;
                int y = (int)mapped_y;
                
                float avg_distance = CALC_AVG_DISTANCE(node);
                
                for(int dy = -2; dy <= 2; dy++) {
                    for(int dx = -2; dx <= 2; dx++) {
                        int px = x + dx;
                        int py = y + dy;
                        
                        if (px >= 0 && px < width && py >= 0 && py < height) {
                            int idx = (py * width + px) * 3;
                            
                            int dist_step = (int)(avg_distance);
                            if (dist_step >= 15) dist_step = 14;
                            
                            image[idx] = (unsigned char)(colors[dist_step][0] * 255);
                            image[idx+1] = (unsigned char)(colors[dist_step][1] * 255);
                            image[idx+2] = (unsigned char)(colors[dist_step][2] * 255);
                        }
                    }
                }
            }
        }
    }

    char* tmp_filename = strdup(filename);
    char* pcd_name = basename(tmp_filename);
    char output_path[1024];
    system("mkdir -p ./output");
    snprintf(output_path, sizeof(output_path), "./output/%s.png", pcd_name);
    
    save_png(output_path, image, width, height);
    free(tmp_filename);
    tracked_free(image, width * height * 3);
    
    double png_end = get_current_time();
    
    metrics.png_generation.time_taken = png_end - png_start;
    metrics.png_generation.memory_used = (peak_memory - before_png_memory) / 1024;  // PNG 생성 중 최대 메모리 사용량
    
    metrics.total.time_taken = png_end - start_time;
    metrics.total.memory_used = peak_memory / 1024;  // 전체 과정 중 최대 메모리 사용량
    
    tracked_free_grid(grid);
    return metrics;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("사용법: %s <pcd_file>\n", argv[0]);
        return 1;
    }

    FILE* result_file = fopen("result.txt", "w");
    if (!result_file) {
        printf("결과 파일을 생성할 수 없습니다.\n");
        return 1;
    }

    ProcessingMetrics all_metrics[NUM_ITERATIONS];
    printf("성능 테스트 시작 (%d회 반복)...\n", NUM_ITERATIONS);
    
    fprintf(result_file, "=== 각 라운드별 성능 측정 결과 ===\n");

    for (int i = 0; i < NUM_ITERATIONS; i++) {
        printf("\r반복 %d/%d", i + 1, NUM_ITERATIONS);
        fflush(stdout);
        reset_memory_tracking();
        all_metrics[i] = process_iteration(argv[1]);

        // 각 라운드의 결과를 파일에 기록
        fprintf(result_file, "\n[Round %d]\n", i + 1);
        fprintf(result_file, "그리드 구축: %.2f ms, %ld KB\n", 
                all_metrics[i].grid_construction.time_taken * 1000.0,
                all_metrics[i].grid_construction.memory_used);
        fprintf(result_file, "PNG 생성: %.2f ms, %ld KB\n",
                all_metrics[i].png_generation.time_taken * 1000.0,
                all_metrics[i].png_generation.memory_used);
        fprintf(result_file, "전체 처리: %.2f ms, %ld KB\n",
                all_metrics[i].total.time_taken * 1000.0,
                all_metrics[i].total.memory_used);
    }
    printf("\n\n");
    
    // 각 단계별 메트릭 수집
    PerfMetrics grid_metrics[NUM_ITERATIONS];
    PerfMetrics png_metrics[NUM_ITERATIONS];
    PerfMetrics total_metrics[NUM_ITERATIONS];
    
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        grid_metrics[i] = all_metrics[i].grid_construction;
        png_metrics[i] = all_metrics[i].png_generation;
        total_metrics[i] = all_metrics[i].total;
    }
    
    // 최종 통계를 파일과 터미널에 출력
    fprintf(result_file, "\n=== 최종 통계 ===\n");
    print_statistics("그리드 구축 단계", grid_metrics, NUM_ITERATIONS, result_file, 1);
    print_statistics("PNG 생성 단계", png_metrics, NUM_ITERATIONS, result_file, 1);
    print_statistics("전체 처리", total_metrics, NUM_ITERATIONS, result_file, 1);
    
    fclose(result_file);
    return 0;
}