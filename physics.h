// シミュレーション条件を格納する構造体
// 反発係数CORを追加
typedef struct condition {
    int width;            // 見えている範囲の幅
    int height;           // 見えている範囲の高さ
    double G;             // 重力定数
    double dt;            // シミュレーションの時間幅
    double cor;           // 壁の反発係数
    double fusion_thresh; // 物体が融合するかどうかの閾値となる2乗距離
} Condition;

// 個々の物体を表す構造体
typedef struct object {
    int exist; // 存在しているか
    double m;
    double x;
    double prev_x;
    double y;
    double prev_y; // 壁からの反発に使用
    double vx;
    double prev_vx; // 壁からの反発に使用
    double vy;
    double prev_vy;
    char symbol[256];
    int nxt_symbol_idx;
} Object;

void my_plot_objects(Object objs[], const size_t numobj, const double t, const Condition cond, const int real_numobj);
void my_update_velocities(Object objs[], const size_t numobj, const Condition cond);
void my_update_positions(Object objs[], const size_t numobj, const Condition cond);
int my_bounce(Object objs[], const size_t numobj, const Condition cond);

int read_data(const char *filename, Object objects[], const size_t numobj, int *real_numobj);
int within_grid(const int y, const int x, const Condition *cond);
void calc_accel(const Object objs[], const size_t numobj, const Condition cond, const int i, double *a_x, double *a_y);
void calc_prev_accel(const Object objs[], const size_t numobj, const Condition cond, const int i, double *a_x,
                     double *a_y);
int collision(Object *obj, const Condition *cond, const double a_x, const double a_y);
void fusion(Object *obj1, Object *obj2);