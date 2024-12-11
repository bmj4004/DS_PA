#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// 설정 가능한 파라미터
#define HORIZONTAL_FOV 70.0f
#define VERTICAL_FOV 43.0f
#define ANGULAR_RESOLUTION 0.5f
#define CENTER_ELEVATION 90.0f
#define CAMERA_HEIGHT 1.0f     // 카메라 높이 (미터 단위)
#define FOV_DISTANCE 3.0f      // FOV 시각화 거리

// OpenGL 관련 전역 변수
float rotation_x = 30.0f;
float rotation_y = -45.0f;
float rotation_z = 0.0f;
float zoom = -20.0f;

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

// 전역 변수로 그리드 선언
SphericalGrid* global_grid = NULL;

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
    HeapNode* node = (HeapNode*)malloc(sizeof(HeapNode));
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

    // 왼쪽, 오른쪽 중 작은 쪽으로 진행
    if (!root->left || (root->right && 
        root->left->distance > root->right->distance)) {
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

// 노드 추가 함수
void add_point_to_grid(SphericalGrid* grid, float x, float y, float z) {
   // 카메라 위치 기준으로 상대 좌표 계산
   float rel_z = z - CAMERA_HEIGHT;
   
   HeapNode* node = create_node(x, y, z);  // 원본 좌표는 보존
   if (!node) return;
   
   // 카메라 기준 각도 계산
   float distance = sqrt(x*x + y*y + rel_z*rel_z);
   float azimuth = atan2(y, x) * 180.0f / M_PI;
   if(azimuth < 0) azimuth += 360.0f;
   float elevation = acos(rel_z / distance) * 180.0f / M_PI;
   
   if (azimuth < 0 || azimuth >= 360.0f ||
       elevation < 0 || elevation >= 180.0f) {
       free(node);
       return;
   }
   
   int azimuth_idx = (int)(azimuth / ANGULAR_RESOLUTION);
   int elevation_idx = (int)(elevation / ANGULAR_RESOLUTION);
   
   if (azimuth_idx < 0 || azimuth_idx >= grid->azimuth_cells ||
       elevation_idx < 0 || elevation_idx >= grid->elevation_cells) {
       free(node);
       return;
   }

   float half_hfov = HORIZONTAL_FOV / 2.0f;
   float half_vfov = VERTICAL_FOV / 2.0f;

   if ((azimuth <= half_hfov || azimuth >= (360.0f - half_hfov)) &&
       (fabs(elevation - CENTER_ELEVATION) <= half_vfov)) {
       
       // 디버그 출력
       printf("Adding point - az: %.2f, el: %.2f [%d][%d] (x:%.2f, y:%.2f, z:%.2f)\n",
              azimuth, elevation, azimuth_idx, elevation_idx, x, y, z);
       
       if (grid->cells[azimuth_idx][elevation_idx].is_valid) {
           node->distance = distance;  // 계산된 거리 저장
           node->azimuth = azimuth;   // 계산된 방위각 저장
           node->elevation = elevation; // 계산된 고도각 저장
           insert_node(&grid->cells[azimuth_idx][elevation_idx].heap, node);
       } else {
           free(node);
       }
   } else {
       free(node);
   }
}

// 구면 그리드 생성
SphericalGrid* create_grid() {
    SphericalGrid* grid = (SphericalGrid*)malloc(sizeof(SphericalGrid));
    if (!grid) return NULL;
    
    grid->azimuth_cells = (int)(360.0f / ANGULAR_RESOLUTION);
    grid->elevation_cells = (int)(180.0f / ANGULAR_RESOLUTION);
    
    grid->cells = (SphericalCell**)malloc(grid->azimuth_cells * sizeof(SphericalCell*));
    if (!grid->cells) {
        free(grid);
        return NULL;
    }

    for (int i = 0; i < grid->azimuth_cells; i++) {
        grid->cells[i] = (SphericalCell*)malloc(grid->elevation_cells * sizeof(SphericalCell));
        if (!grid->cells[i]) {
            for (int j = 0; j < i; j++) {
                free(grid->cells[j]);
            }
            free(grid->cells);
            free(grid);
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
    
    printf("그리드 생성 완료 (크기: %d x %d)\n", grid->azimuth_cells, grid->elevation_cells);
    return grid;
}

// 재귀적으로 힙 노드들을 해제
void free_heap_nodes(HeapNode* node) {
    if (node) {
        free_heap_nodes(node->left);
        free_heap_nodes(node->right);
        free(node);
    }
}

// 그리드 메모리 해제
void free_grid(SphericalGrid* grid) {
    if (!grid) return;

    for (int i = 0; i < grid->azimuth_cells; i++) {
        if (grid->cells[i]) {
            for (int j = 0; j < grid->elevation_cells; j++) {
                free_heap_nodes(grid->cells[i][j].heap.root);
            }
            free(grid->cells[i]);
        }
    }
    free(grid->cells);
    free(grid);
}

// 재귀적으로 힙을 순회하며 그리기
void draw_heap_points(HeapNode* node, float max_distance) {
    if (!node) return;

    float normalized = node->distance / max_distance;
    if (normalized < 0.5f) {
        glColor3f(0.0f, 2.0f * normalized, 1.0f - 2.0f * normalized);
    } else {
        glColor3f(2.0f * (normalized - 0.5f), 2.0f * (1.0f - normalized), 0.0f);
    }
    glVertex3f(node->x, node->y, node->z);

    draw_heap_points(node->left, max_distance);
    draw_heap_points(node->right, max_distance);
}

// 최대 거리 찾기 (재귀)
void find_max_distance(HeapNode* node, float* max_distance) {
    if (!node) return;
    
    if (node->distance > *max_distance) {
        *max_distance = node->distance;
    }
    
    find_max_distance(node->left, max_distance);
    find_max_distance(node->right, max_distance);
}

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

    printf("총 포인트 수: %d\n", num_points);
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
        
        if (i % 10000 == 0) {
            printf("처리 중... %d/%d (유효 포인트: %d)\n", i, num_points, valid_points);
        }
    }
    
    printf("처리 완료. 총 유효 포인트: %d\n", valid_points);
    fclose(file);
    return grid;
}

// OpenGL 관련 함수들
void init(void) {
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(2.0f);
}

void draw_fov() {
    float h_fov_rad = HORIZONTAL_FOV * M_PI / 180.0f;
    float v_fov_rad = VERTICAL_FOV * M_PI / 180.0f;
    
    glColor4f(0.5f, 0.0f, 0.5f, 1.0f);
    glLineWidth(2.0f);
    
    glBegin(GL_LINES);
    
    // 카메라 위치 (높이만큼 올라간 새로운 원점)
    float x0 = 0, y0 = 0, z0 = CAMERA_HEIGHT;
    
    // 왼쪽 아래
    float h1 = -h_fov_rad/2;
    float v1 = -v_fov_rad/2;
    float x1 = FOV_DISTANCE * cos(v1) * cos(h1);
    float y1 = FOV_DISTANCE * cos(v1) * sin(h1);
    float z1 = FOV_DISTANCE * sin(v1) + CAMERA_HEIGHT;
    
    // 왼쪽 위
    float h2 = -h_fov_rad/2;
    float v2 = v_fov_rad/2;
    float x2 = FOV_DISTANCE * cos(v2) * cos(h2);
    float y2 = FOV_DISTANCE * cos(v2) * sin(h2);
    float z2 = FOV_DISTANCE * sin(v2) + CAMERA_HEIGHT;
    
    // 오른쪽 위
    float h3 = h_fov_rad/2;
    float v3 = v_fov_rad/2;
    float x3 = FOV_DISTANCE * cos(v3) * cos(h3);
    float y3 = FOV_DISTANCE * cos(v3) * sin(h3);
    float z3 = FOV_DISTANCE * sin(v3) + CAMERA_HEIGHT;
    
    // 오른쪽 아래
    float h4 = h_fov_rad/2;
    float v4 = -v_fov_rad/2;
    float x4 = FOV_DISTANCE * cos(v4) * cos(h4);
    float y4 = FOV_DISTANCE * cos(v4) * sin(h4);
    float z4 = FOV_DISTANCE * sin(v4) + CAMERA_HEIGHT;
    
    // 원점에서 각 모서리로 이어지는 선
    glVertex3f(x0, y0, z0);
    glVertex3f(x1, y1, z1);
    
    glVertex3f(x0, y0, z0);
    glVertex3f(x2, y2, z2);
    
    glVertex3f(x0, y0, z0);
    glVertex3f(x3, y3, z3);
    
    glVertex3f(x0, y0, z0);
    glVertex3f(x4, y4, z4);
    
    // FOV 끝단을 연결하는 선
    glVertex3f(x1, y1, z1);
    glVertex3f(x2, y2, z2);
    
    glVertex3f(x2, y2, z2);
    glVertex3f(x3, y3, z3);
    
    glVertex3f(x3, y3, z3);
    glVertex3f(x4, y4, z4);
    
    glVertex3f(x4, y4, z4);
    glVertex3f(x1, y1, z1);
    
    glEnd();
}


void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // 시점 변환 - 카메라 높이만큼 시점을 위로 올림
    glTranslatef(0.0f, 0.0f, zoom);
    glRotatef(rotation_x, 1.0f, 0.0f, 0.0f);
    glRotatef(rotation_y, 0.0f, 1.0f, 0.0f);
    glRotatef(rotation_z, 0.0f, 0.0f, 1.0f);
    glTranslatef(0.0f, 0.0f, -CAMERA_HEIGHT);  // 시점을 카메라 높이만큼 위로 이동

    // FOV 시각화
    draw_fov();

    // 좌표축 그리기 - 원점에서 시작
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);  // X축 (빨강)
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(5.0f, 0.0f, 0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);  // Y축 (초록)
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 5.0f, 0.0f);
    glColor3f(0.0f, 0.0f, 1.0f);  // Z축 (파랑)
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 5.0f);
    glEnd();

    // 포인트 클라우드 그리기
    if (global_grid) {
        float max_distance = 0.0f;
        
        // 최대 거리 찾기
        for (int i = 0; i < global_grid->azimuth_cells; i++) {
            for (int j = 0; j < global_grid->elevation_cells; j++) {
                if (global_grid->cells[i][j].is_valid) {
                    find_max_distance(global_grid->cells[i][j].heap.root, &max_distance);
                }
            }
        }

        // 포인트 그리기
        glBegin(GL_POINTS);
        for (int i = 0; i < global_grid->azimuth_cells; i++) {
            for (int j = 0; j < global_grid->elevation_cells; j++) {
                if (global_grid->cells[i][j].is_valid) {
                    draw_heap_points(global_grid->cells[i][j].heap.root, max_distance);
                }
            }
        }
        glEnd();
    }

    glutSwapBuffers();
}

void reshape(int w, int h) {
    if (h == 0) h = 1;
    float aspect = (float)w / (float)h;

    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 27: // ESC
            if (global_grid) {
                free_grid(global_grid);
                global_grid = NULL;
            }
            exit(0);
            break;
        case 'w': rotation_x += 5.0f; break;
        case 's': rotation_x -= 5.0f; break;
        case 'a': rotation_y -= 5.0f; break;
        case 'd': rotation_y += 5.0f; break;
        case 'z': rotation_z += 5.0f; break;
        case 'x': rotation_z -= 5.0f; break;
        case 'q': zoom -= 0.5f; break;
        case 'e': zoom += 0.5f; break;
        case 'r': // 뷰 리셋
            rotation_x = 30.0f;
            rotation_y = -45.0f;
            rotation_z = 0.0f;
            zoom = -20.0f;
            break;
    }
    glutPostRedisplay();
}

void special(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:    rotation_x += 5.0f; break;
        case GLUT_KEY_DOWN:  rotation_x -= 5.0f; break;
        case GLUT_KEY_LEFT:  rotation_y -= 5.0f; break;
        case GLUT_KEY_RIGHT: rotation_y += 5.0f; break;
    }
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("사용법: %s <pcd_file>\n", argv[0]);
        return 1;
    }

    global_grid = process_pcd_file(argv[1]);
    if (!global_grid) {
        printf("PCD 파일 처리 실패\n");
        return 1;
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Point Cloud Viewer with FOV");

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);

    init();

    printf("\n조작 방법:\n");
    printf("W/S 또는 위/아래 화살표: X축 회전\n");
    printf("A/D 또는 좌/우 화살표: Y축 회전\n");
    printf("Z/X: Z축 회전\n");
    printf("Q/E: 줌 인/아웃\n");
    printf("R: 뷰 리셋\n");
    printf("ESC: 종료\n");

    glutMainLoop();

    return 0;
}