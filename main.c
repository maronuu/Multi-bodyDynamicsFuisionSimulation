#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "physics.h"

enum SOLVER_IDX {
    EULER = 0,
    HEUN = 1,
    MIDDLE = 2,
};

//// GLOBAL VARIABLES ////
pthread_mutex_t key_mutex;
int IS_POSE;
int IS_QUIT;
const double FRAME_RATE_LIST_MICRO_SEC[4] = {100000, 10000, 100, 10};
const double FRAME_RATE_PERCENTAGE_LIST_MICRO_SEC[4] = {(double)1 / 1000000, (double)1 / 10000, (double)1 / 1,
                                                        (double)10 / 1};
int FRAME_RATE_IDX;

int within_grid(const int y, const int x, const Condition *cond) {
    return (0 <= y && y < cond->height && 0 <= x && x < cond->width);
}

void my_plot_objects(Object objs[], const size_t numobj, const double t, const Condition cond, const int real_numobj) {
    // origin
    const int origin_y = cond.height / 2;
    const int origin_x = cond.width / 2;
    // (x, y)  --  (x+1, y)
    //    |     (a,b)      |
    // (x, y+1)--(x+1, y+1)
    // where, x <= a < x+1, y <= b < y+1

    // grid
    char grid[cond.height][cond.width];
    // initialize grid
    for (int y = 0; y < cond.height; ++y) {
        for (int x = 0; x < cond.width; ++x) {
            grid[y][x] = 0;
        }
    }
    // fill object
    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;
        // map (obj.x, obj.y) |-> (grid_x, gird-Y)
        const int grid_x = origin_x + floor(objs[obj_idx].x);
        const int grid_y = origin_y + floor(objs[obj_idx].y);
        if (within_grid(grid_y, grid_x, &cond)) {
            grid[grid_y][grid_x] = objs[obj_idx].symbol[0];
        }
    }

    // print field
    // top line
    printf("+");
    for (int i = 0; i < cond.width; ++i)
        printf("-");
    printf("+\r\n");
    for (int y = 0; y < cond.height; ++y) {
        printf("|");
        for (int x = 0; x < cond.width; ++x) {
            if (grid[y][x] != '\0') {
                printf("%c", grid[y][x]);
            } else {
                printf(" ");
            }
        }
        printf("|\r\n");
    }
    printf("+");
    for (int i = 0; i < cond.width; ++i)
        printf("-");
    printf("+\r\n");
    // objects info
    printf("t = %.3f\r\n", t);
    for (int i = 0; i < real_numobj; ++i) {
        if (!objs[i].exist) {
            // dead
            printf("\x1b[31m");
            printf("%c: objs[%d].x = %.2f, objs[%d].y = %.2f", objs[i].symbol[0], i, objs[i].x, i, objs[i].y);
            printf(" includes {");
            for (int j = 0; j < 256; ++j) {
                if (objs[i].symbol[j] == '\0')
                    break;
                printf(" %c ", objs[i].symbol[j]);
            }
            printf("}\x1b[39m\r\n");
        } else {
            printf("%c: objs[%d].x = %.2f, objs[%d].y = %.2f", objs[i].symbol[0], i, objs[i].x, i, objs[i].y);
            printf(" includes {");
            for (int j = 0; j < 256; ++j) {
                if (objs[i].symbol[j] == '\0')
                    break;
                printf(" %c ", objs[i].symbol[j]);
            }
            printf("}\r\n");
        }
    }
}

void calc_accel(const Object objs[], const size_t numobj, const Condition cond, const int i, double *a_x, double *a_y) {
    double ret_y = 0;
    double ret_x = 0;
    for (int j = 0; j < numobj; ++j) {
        if (j == i || !(objs[j].exist) || !(objs[i].exist))
            continue;
        const double dx = objs[j].x - objs[i].x;
        const double dy = objs[j].y - objs[i].y;
        const double dr = sqrt(dx * dx + dy * dy);
        ret_x += objs[j].m * dx / (dr * dr * dr);
        ret_y += objs[j].m * dy / (dr * dr * dr);
    }
    *a_x = cond.G * ret_x;
    *a_y = cond.G * ret_y;
}

void calc_prev_accel(const Object objs[], const size_t numobj, const Condition cond, const int i, double *a_x,
                     double *a_y) {
    double ret_x = 0;
    double ret_y = 0;
    for (int j = 0; j < numobj; ++j) {
        if (j == i || !(objs[j].exist) || !(objs[i].exist))
            continue;
        const double dx = objs[j].prev_x - objs[i].prev_x;
        const double dy = objs[j].prev_y - objs[i].prev_y;
        const double dr = sqrt(dx * dx + dy * dy);
        ret_x += objs[j].m * dx / (dr * dr * dr);
        ret_y += objs[j].m * dy / (dr * dr * dr);
    }
    *a_x = cond.G * ret_x;
    *a_y = cond.G * ret_y;
}

void my_update_velocities(Object objs[], const size_t numobj, const Condition cond) {
    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;
        double a_x, a_y;
        calc_accel(objs, numobj, cond, obj_idx, &a_x, &a_y);
        // x
        objs[obj_idx].prev_vx = objs[obj_idx].vx;
        objs[obj_idx].prev_vy = objs[obj_idx].vy;
        objs[obj_idx].vx += a_x * cond.dt;
        objs[obj_idx].vy += a_y * cond.dt;
    }
}

void my_update_positions(Object objs[], const size_t numobj, const Condition cond) {
    // CAUTION: velocity is already updated
    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;
        objs[obj_idx].prev_x = objs[obj_idx].x;
        objs[obj_idx].prev_y = objs[obj_idx].y;
        objs[obj_idx].x += objs[obj_idx].prev_vx * cond.dt;
        objs[obj_idx].y += objs[obj_idx].prev_vy * cond.dt;
    }
}

// オイラー法による更新
void my_update_params_euler(Object objs[], const size_t numobj, const Condition cond) {
    my_update_velocities(objs, numobj, cond);
    my_update_positions(objs, numobj, cond);
}

// ホイン法による更新
void my_update_params_heuns(Object objs[], const size_t numobj, const Condition cond) {
    double cur_xs[numobj];
    double cur_ys[numobj];
    double cur_vxs[numobj];
    double cur_vys[numobj];
    double cur_axs[numobj];
    double cur_ays[numobj];

    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;
        // calc accel on starting point
        double accel_x_start, accel_y_start;
        calc_accel(objs, numobj, cond, obj_idx, &accel_x_start, &accel_y_start);
        cur_axs[obj_idx] = accel_x_start;
        cur_ays[obj_idx] = accel_y_start;
        // record current
        cur_xs[obj_idx] = objs[obj_idx].x;
        cur_ys[obj_idx] = objs[obj_idx].y;
        cur_vxs[obj_idx] = objs[obj_idx].vx;
        cur_vys[obj_idx] = objs[obj_idx].vy;
    }
    // まず、オイラー法で更新し、t+1の値(暫定)を得る
    my_update_params_euler(objs, numobj, cond);

    // (t+1, y_t+1)(終端)における値
    // double end_xs[numobj];
    // double end_ys[numobj];
    double end_vxs[numobj];
    double end_vys[numobj];
    double end_axs[numobj];
    double end_ays[numobj];

    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;
        // record end point val
        // end_xs[obj_idx] = objs[obj_idx].x;
        // end_ys[obj_idx] = objs[obj_idx].y;
        end_vxs[obj_idx] = objs[obj_idx].vx;
        end_vys[obj_idx] = objs[obj_idx].vy;
        // calc accel
        double accel_x_end, accel_y_end;
        calc_accel(objs, numobj, cond, obj_idx, &accel_x_end, &accel_y_end);
        end_axs[obj_idx] = accel_x_end;
        end_ays[obj_idx] = accel_y_end;
    }

    // update prev
    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;
        objs[obj_idx].prev_x = cur_xs[obj_idx];
        objs[obj_idx].prev_y = cur_ys[obj_idx];
        objs[obj_idx].prev_vx = cur_vxs[obj_idx];
        objs[obj_idx].prev_vy = cur_vys[obj_idx];
    }

    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;

        // update v
        objs[obj_idx].vx = cur_vxs[obj_idx] + (cond.dt / 2.0) * (cur_axs[obj_idx] + end_axs[obj_idx]);
        objs[obj_idx].vy = cur_vys[obj_idx] + (cond.dt / 2.0) * (cur_ays[obj_idx] + end_ays[obj_idx]);
        // update r
        objs[obj_idx].x = cur_xs[obj_idx] + (cond.dt / 2.0) * (cur_vxs[obj_idx] + end_vxs[obj_idx]);
        objs[obj_idx].y = cur_ys[obj_idx] + (cond.dt / 2.0) * (cur_vys[obj_idx] + end_vys[obj_idx]);
    }
}

// 中点法による更新
void my_update_params_middle_point_method(Object objs[], const size_t numobj, const Condition cond) {
    double cur_xs[numobj];
    double cur_ys[numobj];
    double cur_vxs[numobj];
    double cur_vys[numobj];

    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;
        // record current
        cur_xs[obj_idx] = objs[obj_idx].x;
        cur_ys[obj_idx] = objs[obj_idx].y;
        cur_vxs[obj_idx] = objs[obj_idx].vx;
        cur_vys[obj_idx] = objs[obj_idx].vy;
        // calc r on middle
        const double rx_mid = objs[obj_idx].x + (cond.dt / 2.0) * objs[obj_idx].vx;
        const double ry_mid = objs[obj_idx].y + (cond.dt / 2.0) * objs[obj_idx].vy;
        objs[obj_idx].x = rx_mid;
        objs[obj_idx].y = ry_mid;
    }

    // update prev
    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        objs[obj_idx].prev_x = cur_xs[obj_idx];
        objs[obj_idx].prev_y = cur_ys[obj_idx];
        objs[obj_idx].prev_vx = cur_vxs[obj_idx];
        objs[obj_idx].prev_vy = cur_vys[obj_idx];
    }

    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;

        // calc accel on middle
        double accel_x_mid, accel_y_mid;
        calc_accel(objs, numobj, cond, obj_idx, &accel_x_mid, &accel_y_mid);
        // update v
        objs[obj_idx].vx = objs[obj_idx].vx + accel_x_mid * cond.dt;
        objs[obj_idx].vy = objs[obj_idx].vy + accel_y_mid * cond.dt;
        // update x
        double vx_mid = cur_vxs[obj_idx] + (cond.dt / 2.0) * accel_x_mid;
        double vy_mid = cur_vys[obj_idx] + (cond.dt / 2.0) * accel_y_mid;
        objs[obj_idx].x = cur_xs[obj_idx] + vx_mid * cond.dt;
        objs[obj_idx].y = cur_ys[obj_idx] + vy_mid * cond.dt;
    }
}

int collision(Object *obj, const Condition *cond, const double a_x, const double a_y) {
    const double origin_y = (double)(cond->height / 2);
    const double origin_x = (double)(cond->width / 2);

    const double top = 0.0 - origin_y;
    const double bottom = cond->height - origin_y;
    const double left = 0.0 - origin_x;
    const double right = cond->width - origin_x;

    double wall;
    int is_y_bounce = 0; // y方向の反発なら1, そうでなければ0
    if (obj->y < top && obj->vy <= 0) {
        wall = top;
        is_y_bounce = 1;
    } else if (bottom <= obj->y && obj->vy >= 0) {
        wall = bottom;
        is_y_bounce = 1;
    } else if (obj->x < left && obj->vx <= 0) {
        wall = left;
    } else if (right <= obj->x && obj->vx >= 0) {
        wall = right;
    } else {
        return 0;
    }

    if (is_y_bounce) {
        // before bouncing
        const double time_to_wall = (wall - obj->prev_y) / obj->prev_vy;
        // printf("collide y\n");
        if (!(0 <= time_to_wall && time_to_wall < cond->dt)) {
            fprintf(stderr, "To simulate more extreme dynamics, give more precise dt\r\n");
            system("/bin/stty cooked");
            exit(EXIT_FAILURE);
        }
        const double original_y = obj->y;
        const double original_vy = obj->vy;
        obj->vy = obj->prev_vy + a_y * time_to_wall;
        // bounce! -> correct velocity
        obj->vy *= -cond->cor;
        // after bouncing
        const double remain_time = cond->dt - time_to_wall;
        obj->y = wall + (obj->vy) * remain_time;
        obj->vy += a_y * remain_time;

        obj->prev_y = original_y;
        obj->prev_vy = original_vy;
    } else {
        const double time_to_wall = (wall - obj->prev_x) / obj->prev_vx;
        // printf("collide x\n");
        if (!(0 <= time_to_wall && time_to_wall < cond->dt)) {
            fprintf(stderr, "To simulate more extreme dynamics, give more precise dt\r\n");
            system("/bin/stty cooked");
            exit(EXIT_FAILURE);
        }
        const double original_x = obj->x;
        const double original_vx = obj->vx;
        obj->vx = obj->prev_vx + a_x * cond->dt;
        // bounce! -> correct velocity
        obj->vx *= -cond->cor;
        // after bouncing
        const double remain_time = cond->dt - time_to_wall;
        obj->x = wall + (-obj->prev_vx) * remain_time;

        obj->prev_x = original_x;
        obj->prev_vx = original_vx;
    }
    return 1;
}

int my_bounce(Object objs[], const size_t numobj, const Condition cond) {
    int is_any_bounce = 0;
    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        if (!objs[obj_idx].exist)
            continue;
        // calc acceleration for time = t
        double prev_a_x, prev_a_y;
        calc_prev_accel(objs, numobj, cond, obj_idx, &prev_a_x, &prev_a_y);
        if (collision(&objs[obj_idx], &cond, prev_a_x, prev_a_y)) {
            is_any_bounce = 1;
        }
    }
    return is_any_bounce;
}

double calc_sq_distance(const Object *obj1, const Object *obj2) {
    const double dx = obj1->x - obj2->x;
    const double dy = obj1->y - obj2->y;
    return dx * dx + dy * dy;
}

void fusion(Object *obj1, Object *obj2) {
    const double m1 = obj1->m;
    const double m2 = obj2->m;
    const double vx = (m1 * obj1->vx + m2 * obj2->vx) / (m1 + m2);
    const double vy = (m1 * obj1->vy + m2 * obj2->vy) / (m1 + m2);

    obj1->m = m1 + m2;
    obj1->vx = vx;
    obj1->vy = vy;

    // symbol merge
    obj1->symbol[obj1->nxt_symbol_idx++] = obj2->symbol[0];
    obj2->symbol[1] = '\0';

    obj2->exist = 0;
}

int fusion_pairs(Object objs[], const size_t numobj, const Condition cond) {
    int is_any_fusion = 0;
    for (int i = 0; i < numobj; ++i) {
        if (!objs[i].exist)
            continue;
        for (int j = i + 1; j < numobj; ++j) {
            if (!objs[j].exist)
                continue;
            if (calc_sq_distance(&objs[i], &objs[j]) < cond.fusion_thresh) {
                fusion(&objs[i], &objs[j]);
                is_any_fusion = 1;
            }
        }
    }
    return is_any_fusion;
}

// thread function to wait key input
void *wait_key() {
    char ch;
    while (1) {
        ch = getchar();
        pthread_mutex_lock(&key_mutex);
        switch (ch) {
        case 'q':
            IS_QUIT = 1;
            printf("\r\nSimulation is over\r\n");
            system("/bin/stty cooked");
            exit(0);
        // case 'f':
        //     IS_POSE = 1;
        //     printf("Press 'f' to restart\r\n");
        //     printf("In thread\r\n");
        //     char key;
        //     while ((key = getchar()) != 'f')
        //         ;
        //     IS_POSE = 0;
        //     break;
        case 'd': // speed down
            if (FRAME_RATE_IDX < 3)
                FRAME_RATE_IDX++;
            break;
        case 'a': // speed up
            if (FRAME_RATE_IDX > 0)
                FRAME_RATE_IDX--;
            break;
        default:
            break;
        }
        pthread_mutex_unlock(&key_mutex);
        usleep(100);
    };

    return (void *)ch;
}

int read_data(const char *filename, Object objects[], const size_t numobj, int *real_numobj) {
    // initialize all objs
    for (int obj_idx = 0; obj_idx < numobj; ++obj_idx) {
        objects[obj_idx].exist = 0;
        objects[obj_idx].m = 0;
    }

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Cannot Open file\n");
        return 1;
    }
    char buffer[512];
    while (fgets(buffer, 512, fp)) {
        if (buffer[0] != '#')
            break;
    }
    int obj_idx = 0;
    int used_symbol[128];
    for (int i = 0; i < 128; ++i)
        used_symbol[i] = 0;
    while (buffer != NULL) {
        int i = 0;
        while (buffer[i] != '#' && buffer[i] != '\0' && buffer[i] != '\n')
            ++i;
        buffer[i] = '\0';
        --i;
        while (buffer[i] == ' ') {
            buffer[i] = '\0';
            --i;
        }
        double m, x, y, vx, vy;
        sscanf(buffer, "%lf %lf %lf %lf %lf", &m, &x, &y, &vx, &vy);
        // printf("m = %lf, x = %lf, y = %lf, vx = %lf, vy = %lf\n", m, x, y, vx, vy);
        objects[obj_idx].exist = 1;
        objects[obj_idx].m = m;
        objects[obj_idx].x = x;
        objects[obj_idx].y = y;
        objects[obj_idx].vx = vx;
        objects[obj_idx].vy = vy;
        // choose symbol randomly
        int diff = rand() % 26;
        char sym = 'A' + diff;
        while (used_symbol[diff]) {
            diff = rand() % 26;
            sym = 'A' + diff;
        }
        used_symbol[diff] = 1;
        for (int si = 0; si < 256; ++si)
            objects[obj_idx].symbol[si] = '\0';
        objects[obj_idx].symbol[0] = sym;
        objects[obj_idx].nxt_symbol_idx = 1;

        obj_idx++;
        if (obj_idx >= numobj)
            break;
        if (fgets(buffer, 512, fp) == NULL)
            break;
    }
    *real_numobj = obj_idx;
    fclose(fp);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 5) {
        fprintf(stderr, "Too many args\n");
        fprintf(stderr, "Usage: ./main 5 path/to/datafile.dat middle");
        return EXIT_FAILURE;
    } else if (argc <= 3) {
        fprintf(stderr, "Too few args\n");
        fprintf(stderr, "Usage: ./main 5 path/to/datafile.dat middle");
        return EXIT_FAILURE;
    }
    system("/bin/stty raw onlcr"); // enterを押さなくてもキー入力を受け付けるようになる

    // set seed
    srand(time(NULL));

    const Condition cond = {
        .width = 75,
        .height = 40,
        .G = 1.0,
        .dt = 5e-2,
        .cor = 0.8,
        .fusion_thresh = 2,
    };

    int objnum_ = atoi(argv[1]);
    if (objnum_ < 1) {
        fprintf(stderr, "number of objects should be >=1\r\n");
        return EXIT_FAILURE;
    }
    size_t objnum = (size_t)objnum_;
    Object objects[objnum];

    int real_objnum;
    if (read_data(argv[2], objects, objnum, &real_objnum) == 1) {
        fprintf(stderr, "File Error\r\n");
        return EXIT_FAILURE;
    }

    // solver selection
    int solver_idx;
    if (strcmp(argv[3], "euler") == 0) {
        solver_idx = EULER;
    } else if (strcmp(argv[3], "heun") == 0) {
        solver_idx = HEUN;
    } else if (strcmp(argv[3], "middle") == 0) {
        solver_idx = MIDDLE;
    } else {
        fprintf(stderr, "Invalid Solver %s\r\n", argv[3]);
        fprintf(stderr, "Please select ('euler' or 'heun' or 'middle')");
        return EXIT_FAILURE;
    }

    // thread
    pthread_t wait_key_thread;
    int thread_status = pthread_create(&wait_key_thread, NULL, (void *)wait_key, NULL);
    if (thread_status != 0) {
        fprintf(stderr, "Cannot create thread\r\n");
        return EXIT_FAILURE;
    }
    pthread_mutex_init(&key_mutex, NULL);

    // シミュレーション. ループは整数で回しつつ、実数時間も更新する
    FRAME_RATE_IDX = 2;
    const double stop_time = 400;
    double t = 0;
    for (size_t i = 0; t <= stop_time; i++) {
        system("clear");
        usleep(1000);
        printf("Press [q] to end simulation\r\n");
        printf("Speed: [a] <<");
        for (int f_id = 0; f_id < 4; ++f_id) {
            if (f_id == FRAME_RATE_IDX) {
                printf(" \e[32m[[%fx]]\e[39m ", FRAME_RATE_PERCENTAGE_LIST_MICRO_SEC[f_id]);
            } else {
                printf(" %fx ", FRAME_RATE_PERCENTAGE_LIST_MICRO_SEC[f_id]);
            }
        }
        printf("<< [d]\r\n");

        t = i * cond.dt;

        if (solver_idx == EULER) {
            // オイラー法による更新
            my_update_params_euler(objects, objnum, cond);
        } else if (solver_idx == HEUN) {
            // ホイン法による更新
            my_update_params_heuns(objects, objnum, cond);
        } else if (solver_idx == MIDDLE) {
            // 中点法による更新
            my_update_params_middle_point_method(objects, objnum, cond);
        } else {
            fprintf(stderr, "Invalid Solver idx %d\n", solver_idx);
            system("/bin/stty cooked");
            return EXIT_FAILURE;
        }

        while (my_bounce(objects, objnum, cond) == 1) {
        }
        while (fusion_pairs(objects, objnum, cond) == 1) {
        };

        if (IS_QUIT) {
            break;
        }

        // 表示の座標系は width/2, height/2 のピクセル位置が原点となるようにする
        my_plot_objects(objects, objnum, t, cond, real_objnum);

        // debug
        // printf("time=%.1f, y=%.2f, prev_y=%.2f, vy=%.2f, prev_vy=%.2f\n",
        //  t, objects[0].y, objects[0].prev_y, objects[0].vy, objects[0].prev_vy);

        usleep(FRAME_RATE_LIST_MICRO_SEC[FRAME_RATE_IDX]);
        const int curdiff = 6 + real_objnum;
        printf("\e[%dA", cond.height + curdiff);
        usleep(FRAME_RATE_LIST_MICRO_SEC[FRAME_RATE_IDX]);
    }

    system("clear");
    system("/bin/stty cooked");

    pthread_mutex_destroy(&key_mutex);

    printf("Multibody simulation is over!\r\n");

    return EXIT_SUCCESS;
}
