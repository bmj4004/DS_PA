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
#include <float.h>

// 설정 가능한 파라미터
#define HORIZONTAL_FOV 70.0f
#define VERTICAL_FOV 43.0f
#define CAMERA_HEIGHT 1.0f
#define NUM_ITERATIONS 100
#define VOXEL_SIZE 0.1f  // 10cm 단위의 복셀

// 성능 측정을 위한 구조체
typedef struct {
    double time_taken;
    long memory_used;
} PerfMetrics;

typedef struct {
    PerfMetrics grid_construction;
    PerfMetrics png_generation;
    PerfMetrics total;
} ProcessingMetrics;

// 포인트 구조체
typedef struct {
    float x, y, z;
    float distance;  // 원점으로부터의 거리
} Point;

// 3D 그리드 구조체
typedef struct {
    Point*** data;      // 3D 배열
    int size_x;         // x 방향 크기
    int size_y;         // y 방향 크기
    int size_z;         // z 방향 크기
    float min_x;        // x 최소값
    float min_y;        // y 최소값
    float min_z;        // z 최소값
    float max_x;        // x 최대값
    float max_y;        // y 최대값
    float max_z;        // z 최대값
} SpatialGrid;

// 전역 변수로 메모리 추적
static long current_memory = 0;
static long peak_memory = 0;

// 색상 배열 정의
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

// 메모리 관리 함수들
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

void tracked_free(void* ptr, size_t size) {
    if (ptr) {
        current_memory -= size;
        free(ptr);
    }
}

long get_memory_usage() {
    return current_memory / 1024;  // bytes to KB
}

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

// PCD 파일에서 범위 찾기
void find_pcd_bounds(const char* filename, float* min_x, float* min_y, float* min_z,
                    float* max_x, float* max_y, float* max_z) {
    FILE* file = fopen(filename, "rb");
    if (!file) return;

    *min_x = *min_y = *min_z = FLT_MAX;
    *max_x = *max_y = *max_z = -FLT_MAX;

    char line[1024];
    int num_points = 0;
    // 헤더 파싱
    while (fgets(line, sizeof(line), file)) {
        if (strncmp(line, "POINTS", 6) == 0) {
            sscanf(line, "POINTS %d", &num_points);
        }
        else if (strcmp(line, "DATA binary\n") == 0) {
            break;
        }
    }

    for (int i = 0; i < num_points; i++) {
        float x = read_float_le(file);
        float y = read_float_le(file);
        float z = read_float_le(file);
        unsigned char intensity;
        fread(&intensity, 1, 1, file);

        if (fabs(x) < 100.0f && fabs(y) < 100.0f && fabs(z) < 100.0f) {
            *min_x = fminf(*min_x, x);
            *min_y = fminf(*min_y, y);
            *min_z = fminf(*min_z, z);
            *max_x = fmaxf(*max_x, x);
            *max_y = fmaxf(*max_y, y);
            *max_z = fmaxf(*max_z, z);
        }
    }

    fclose(file);
}

// 3D 그리드 생성
SpatialGrid* create_spatial_grid(float min_x, float min_y, float min_z,
                               float max_x, float max_y, float max_z) {
    SpatialGrid* grid = (SpatialGrid*)tracked_malloc(sizeof(SpatialGrid));
    if (!grid) return NULL;

    // 범위를 VOXEL_SIZE로 나누어 배열 크기 계산
    grid->size_x = (int)ceil((max_x - min_x) / VOXEL_SIZE) + 1;
    grid->size_y = (int)ceil((max_y - min_y) / VOXEL_SIZE) + 1;
    grid->size_z = (int)ceil((max_z - min_z) / VOXEL_SIZE) + 1;

    grid->min_x = min_x;
    grid->min_y = min_y;
    grid->min_z = min_z;
    grid->max_x = max_x;
    grid->max_y = max_y;
    grid->max_z = max_z;

    // 3D 배열 할당
    grid->data = (Point***)tracked_malloc(grid->size_x * sizeof(Point**));
    if (!grid->data) {
        tracked_free(grid, sizeof(SpatialGrid));
        return NULL;
    }

    for (int i = 0; i < grid->size_x; i++) {
        grid->data[i] = (Point**)tracked_malloc(grid->size_y * sizeof(Point*));
        if (!grid->data[i]) {
            for (int j = 0; j < i; j++) {
                tracked_free(grid->data[j], grid->size_y * sizeof(Point*));
            }
            tracked_free(grid->data, grid->size_x * sizeof(Point**));
            tracked_free(grid, sizeof(SpatialGrid));
            return NULL;
        }

        for (int j = 0; j < grid->size_y; j++) {
            grid->data[i][j] = (Point*)tracked_malloc(grid->size_z * sizeof(Point));
            if (!grid->data[i][j]) {
                for (int k = 0; k < j; k++) {
                    tracked_free(grid->data[i][k], grid->size_z * sizeof(Point));
                }
                for (int k = 0; k < i; k++) {
                    for (int l = 0; l < grid->size_y; l++) {
                        tracked_free(grid->data[k][l], grid->size_z * sizeof(Point));
                    }
                    tracked_free(grid->data[k], grid->size_y * sizeof(Point*));
                }
                tracked_free(grid->data[i], grid->size_y * sizeof(Point*));
                tracked_free(grid->data, grid->size_x * sizeof(Point**));
                tracked_free(grid, sizeof(SpatialGrid));
                return NULL;
            }
            
            // 초기화
            for (int k = 0; k < grid->size_z; k++) {
                grid->data[i][j][k].distance = FLT_MAX;  // 아직 포인트가 없음을 표시
            }
        }
    }

    return grid;
}

// 그리드에 포인트 추가
void add_point_to_grid(SpatialGrid* grid, float x, float y, float z) {
    int idx_x = (int)((x - grid->min_x) / VOXEL_SIZE);
    int idx_y = (int)((y - grid->min_y) / VOXEL_SIZE);
    int idx_z = (int)((z - grid->min_z) / VOXEL_SIZE);

    if (idx_x < 0 || idx_x >= grid->size_x ||
        idx_y < 0 || idx_y >= grid->size_y ||
        idx_z < 0 || idx_z >= grid->size_z) {
        return;
    }

    float distance = sqrt(x*x + y*y + (z-CAMERA_HEIGHT)*(z-CAMERA_HEIGHT));

    // 이미 있는 포인트보다 더 가까우면 업데이트
    if (distance < grid->data[idx_x][idx_y][idx_z].distance) {
        grid->data[idx_x][idx_y][idx_z].x = x;
        grid->data[idx_x][idx_y][idx_z].y = y;
        grid->data[idx_x][idx_y][idx_z].z = z;
        grid->data[idx_x][idx_y][idx_z].distance = distance;
    }
}

// 그리드 메모리 해제
void free_spatial_grid(SpatialGrid* grid) {
    if (!grid) return;

    for (int i = 0; i < grid->size_x; i++) {
        for (int j = 0; j < grid->size_y; j++) {
            tracked_free(grid->data[i][j], grid->size_z * sizeof(Point));
        }
        tracked_free(grid->data[i], grid->size_y * sizeof(Point*));
    }
    tracked_free(grid->data, grid->size_x * sizeof(Point**));
    tracked_free(grid, sizeof(SpatialGrid));
}

// PCD 파일 처리
SpatialGrid* process_pcd_file(const char* filename) {
    float min_x, min_y, min_z, max_x, max_y, max_z;
    find_pcd_bounds(filename, &min_x, &min_y, &min_z, &max_x, &max_y, &max_z);

    SpatialGrid* grid = create_spatial_grid(min_x, min_y, min_z, max_x, max_y, max_z);
    if (!grid) return NULL;

    FILE* file = fopen(filename, "rb");
    if (!file) {
        free_spatial_grid(grid);
        return NULL;
    }

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

    for (int i = 0; i < num_points; i++) {
        float x = read_float_le(file);
        float y = read_float_le(file);
        float z = read_float_le(file);
        unsigned char intensity;
        fread(&intensity, 1, 1, file);
        
        if (fabs(x) < 100.0f && fabs(y) < 100.0f && fabs(z) < 100.0f) {
            add_point_to_grid(grid, x, y, z);
        }
    }

    fclose(file);
    return grid;
}

// Range view 이미지 생성
void generate_range_view(SpatialGrid* grid, const char* output_path) {
    const int width = 1920;
    const int height = 1080;
    unsigned char* image = (unsigned char*)tracked_malloc(width * height * 3);
    memset(image, 0, width * height * 3);

    float half_hfov = HORIZONTAL_FOV / 2.0f;
    float half_vfov = VERTICAL_FOV / 2.0f;

    // 모든 포인트에 대해
    for (int i = 0; i < grid->size_x; i++) {
        for (int j = 0; j < grid->size_y; j++) {
            for (int k = 0; k < grid->size_z; k++) {
                if (grid->data[i][j][k].distance != FLT_MAX) {
                    Point p = grid->data[i][j][k];
                    
                    // 방위각과 고도각 계산
                    float azimuth = atan2(p.y, p.x) * 180.0f / M_PI;
                    if (azimuth < 0) azimuth += 360.0f;
                    
                    float elevation = acos((p.z - CAMERA_HEIGHT) / p.distance) * 180.0f / M_PI;

                    // FOV 내부인지 확인
                    if ((azimuth <= half_hfov || azimuth >= (360.0f - half_hfov)) &&
                        (fabs(elevation - 90.0f) <= half_vfov)) {
                        
                        // 이미지 좌표로 변환
                        float mapped_x;
                        if (azimuth >= 324.0f) {
                            mapped_x = (azimuth - 324.0f) / (36.0f) * (width/2);
                        } else {
                            mapped_x = (width/2) + (azimuth / 36.0f) * (width/2);
                        }
                        
                        float mapped_y = ((111.0f - elevation) / (111.0f - 69.0f)) * height;
                        
                        int x = width - 1 - (int)mapped_x;
                        int y = (int)mapped_y;

                        // 점 크기 확장하여 그리기
                        for(int dy = -2; dy <= 2; dy++) {
                            for(int dx = -2; dx <= 2; dx++) {
                                int px = x + dx;
                                int py = y + dy;
                                
                                if (px >= 0 && px < width && py >= 0 && py < height) {
                                    int idx = (py * width + px) * 3;
                                    
                                    int dist_step = (int)(p.distance);
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
        }
    }

    save_png(output_path, image, width, height);
    tracked_free(image, width * height * 3);
}

// 단일 반복 처리
ProcessingMetrics process_iteration(const char* filename) {
    ProcessingMetrics metrics = {0};
    reset_memory_tracking();
    
    double start_time = get_current_time();
    
    // 그리드 구축 단계
    double grid_start = get_current_time();
    SpatialGrid* grid = process_pcd_file(filename);
    double grid_end = get_current_time();
    
    metrics.grid_construction.time_taken = grid_end - grid_start;
    metrics.grid_construction.memory_used = get_memory_usage();

    if (!grid) {
        printf("PCD 파일 처리 실패\n");
        return metrics;
    }
    
    // PNG 생성 단계
    double png_start = get_current_time();
    long before_png_memory = current_memory;

    char* tmp_filename = strdup(filename);
    char* pcd_name = basename(tmp_filename);
    char output_path[1024];
    system("mkdir -p ./output");
    snprintf(output_path, sizeof(output_path), "./output/%s.png", pcd_name);
    
    generate_range_view(grid, output_path);
    free(tmp_filename);
    
    double png_end = get_current_time();
    
    metrics.png_generation.time_taken = png_end - png_start;
    metrics.png_generation.memory_used = (peak_memory - before_png_memory) / 1024;
    
    metrics.total.time_taken = png_end - start_time;
    metrics.total.memory_used = peak_memory / 1024;
    
    free_spatial_grid(grid);
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