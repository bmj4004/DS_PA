#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    float x;
    float y;
    float z;
    unsigned char intensity;
} Point;

// 전역 변수
Point* points = NULL;
int num_points = 0;
float rotation_x = 30.0f;
float rotation_y = -45.0f;
float rotation_z = 0.0f;  // z축 회전 추가
float zoom = -20.0f;

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

// PCD 파일 읽기 함수
Point* read_pcd_file(const char* filename, int* point_count) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        printf("파일을 열 수 없습니다: %s\n", filename);
        return NULL;
    }

    // 헤더 파싱
    char line[1024];
    int num_points = 0;
    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;
        
        if (strncmp(line, "POINTS", 6) == 0) {
            sscanf(line, "POINTS %d", &num_points);
        }
        else if (strcmp(line, "DATA binary") == 0) {
            break;
        }
    }

    if (num_points == 0) {
        printf("포인트 수를 읽을 수 없습니다.\n");
        fclose(file);
        return NULL;
    }

    // 메모리 할당
    Point* points = (Point*)malloc(num_points * sizeof(Point));
    if (!points) {
        printf("메모리 할당 실패\n");
        fclose(file);
        return NULL;
    }

    // 데이터 범위 추적 변수
    float min_x = 1e10f, max_x = -1e10f;
    float min_y = 1e10f, max_y = -1e10f;
    float min_z = 1e10f, max_z = -1e10f;

    // 데이터 읽기
    int valid_count = 0;
    for (int i = 0; i < num_points; i++) {
        float x = read_float_le(file);
        float y = read_float_le(file);
        float z = read_float_le(file);
        unsigned char intensity;
        fread(&intensity, 1, 1, file);

        // 유효한 데이터만 저장
        if ((fabs(x) > 0.001f || fabs(y) > 0.001f || fabs(z) > 0.001f) &&
            fabs(x) < 100.0f && fabs(y) < 100.0f && fabs(z) < 100.0f) {
            
            points[valid_count].x = x;
            points[valid_count].y = y;
            points[valid_count].z = z;
            points[valid_count].intensity = intensity;

            // 범위 업데이트
            min_x = fmin(min_x, x); max_x = fmax(max_x, x);
            min_y = fmin(min_y, y); max_y = fmax(max_y, y);
            min_z = fmin(min_z, z); max_z = fmax(max_z, z);

            valid_count++;
        }

        if (i < 5 || i % 10000 == 0) {
            printf("Point %d: x=%.3f, y=%.3f, z=%.3f, intensity=%d\n",
                   i, x, y, z, intensity);
        }
    }

    printf("\n데이터 범위:\n");
    printf("X: %.3f ~ %.3f\n", min_x, max_x);
    printf("Y: %.3f ~ %.3f\n", min_y, max_y);
    printf("Z: %.3f ~ %.3f\n", min_z, max_z);

    *point_count = valid_count;
    fclose(file);
    return points;
}

// 높이에 따른 색상 계산 함수
void get_height_color(float z, float min_z, float max_z, float* r, float* g, float* b) {
    float normalized_z = (z - min_z) / (max_z - min_z);
    if (normalized_z < 0.5f) {
        *r = 0.0f;
        *g = 2.0f * normalized_z;
        *b = 1.0f - 2.0f * normalized_z;
    } else {
        *r = 2.0f * (normalized_z - 0.5f);
        *g = 2.0f * (1.0f - normalized_z);
        *b = 0.0f;
    }
}

void init(void) {
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(2.0f);
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // 시점 변환
    glTranslatef(0.0f, 0.0f, zoom);
    glRotatef(rotation_x, 1.0f, 0.0f, 0.0f);
    glRotatef(rotation_y, 0.0f, 1.0f, 0.0f);
    glRotatef(rotation_z, 0.0f, 0.0f, 1.0f);  // z축 회전 추가

    // 좌표축 그리기
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f); // X축 (빨강)
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(5.0f, 0.0f, 0.0f);
    glColor3f(0.0f, 1.0f, 0.0f); // Y축 (초록)
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 5.0f, 0.0f);
    glColor3f(0.0f, 0.0f, 1.0f); // Z축 (파랑)
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 5.0f);
    glEnd();

    // 포인트 클라우드 그리기
    glBegin(GL_POINTS);
    
    // z 범위 계산
    float min_z = 1e10f, max_z = -1e10f;
    for (int i = 0; i < num_points; i++) {
        if (points[i].z < min_z) min_z = points[i].z;
        if (points[i].z > max_z) max_z = points[i].z;
    }

    for (int i = 0; i < num_points; i++) {
        float r, g, b;
        get_height_color(points[i].z, min_z, max_z, &r, &g, &b);
        glColor3f(r, g, b);
        glVertex3f(points[i].x, points[i].y, points[i].z);
    }
    glEnd();

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
            if (points) free(points);
            exit(0);
            break;
        case 'w': rotation_x += 5.0f; break;
        case 's': rotation_x -= 5.0f; break;
        case 'a': rotation_y -= 5.0f; break;
        case 'd': rotation_y += 5.0f; break;
        case 'z': rotation_z += 5.0f; break;  // z축 시계 방향
        case 'x': rotation_z -= 5.0f; break;  // z축 반시계 방향
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

    points = read_pcd_file(argv[1], &num_points);
    if (!points) {
        printf("PCD 파일 로드 실패\n");
        return 1;
    }

    printf("로드된 포인트 수: %d\n", num_points);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Point Cloud Viewer");

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