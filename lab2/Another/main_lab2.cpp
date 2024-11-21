#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <locale> 
 

const double EPS = 1e-6;
const int MAX_ITER = 1000;

// Функция для загрузки матрицы в диагональном формате
void read_matrix_diagonal_format(const std::string& filename, std::vector<double>& Dm2, std::vector<double>& Dm1,
                                 std::vector<double>& D0, std::vector<double>& D1, std::vector<double>& D2, int& n) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла матрицы!\n";
        return;
    }
    file >> n;
    Dm2.resize(n - 2);
    Dm1.resize(n - 1);
    D0.resize(n);
    D1.resize(n - 1);
    D2.resize(n - 2);
    for (int i = 0; i < n - 2; ++i) file >> Dm2[i];
    for (int i = 0; i < n - 1; ++i) file >> Dm1[i];
    for (int i = 0; i < n; ++i) file >> D0[i];
    for (int i = 0; i < n - 1; ++i) file >> D1[i];
    for (int i = 0; i < n - 2; ++i) file >> D2[i];
}

// Функция для загрузки правой части и начального приближения
void read_vector(const std::string& filename, std::vector<double>& vec, int n) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла!\n";
        return;
    }
    vec.resize(n);
    for (int i = 0; i < n; ++i) {
        file >> vec[i];
    }
}

// Подсчет невязки
double calculate_residual(const std::vector<double>& x, const std::vector<double>& b, const std::vector<double>& Dm2,
                          const std::vector<double>& Dm1, const std::vector<double>& D0, const std::vector<double>& D1,
                          const std::vector<double>& D2, int n) {
    double norm_b = 0.0, norm_r = 0.0;
    for (int i = 0; i < n; ++i) {
        double Ax_i = D0[i] * x[i];
        if (i > 0) Ax_i += Dm1[i - 1] * x[i - 1];
        if (i < n - 1) Ax_i += D1[i] * x[i + 1];
        if (i > 1) Ax_i += Dm2[i - 2] * x[i - 2];
        if (i < n - 2) Ax_i += D2[i] * x[i + 2];
        double r_i = b[i] - Ax_i;
        norm_b += b[i] * b[i];
        norm_r += r_i * r_i;
    }
    return std::sqrt(norm_r / norm_b);
}

// Итерационный шаг для методов Якоби и Гаусса-Зейделя
void iteration_step(std::vector<double>& x, const std::vector<double>& b, const std::vector<double>& Dm2,
                    const std::vector<double>& Dm1, const std::vector<double>& D0, const std::vector<double>& D1,
                    const std::vector<double>& D2, int n, double omega, bool gauss_seidel) {
    std::vector<double> new_x = x;
    for (int i = 0; i < n; ++i) {
        double sum = b[i];
        if (i > 0) sum -= Dm1[i - 1] * (gauss_seidel ? new_x[i - 1] : x[i - 1]);
        if (i < n - 1) sum -= D1[i] * x[i + 1];
        if (i > 1) sum -= Dm2[i - 2] * (gauss_seidel ? new_x[i - 2] : x[i - 2]);
        if (i < n - 2) sum -= D2[i] * x[i + 2];
        new_x[i] = (1 - omega) * x[i] + omega * (sum / D0[i]);
    }
    x = new_x;
}

// Основная функция решения СЛАУ
void solve_SLAE(const std::string& config_file, const std::string& matrix_file, const std::string& rhs_file,
                const std::string& initial_guess_file, const std::string& output_file, bool gauss_seidel, double omega) {
    int n;
    std::vector<double> Dm2, Dm1, D0, D1, D2;
    read_matrix_diagonal_format(matrix_file, Dm2, Dm1, D0, D1, D2, n);

    std::vector<double> b, x;
    read_vector(rhs_file, b, n);
    read_vector(initial_guess_file, x, n);

    // Итерационный процесс
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        iteration_step(x, b, Dm2, Dm1, D0, D1, D2, n, omega, gauss_seidel);

        double residual = calculate_residual(x, b, Dm2, Dm1, D0, D1, D2, n);
        std::cout << "Iteration " << iter + 1 << ": residual = " << residual << "\n";

        if (residual < EPS) {
            std::cout << "Сходимость достигнута на итерации " << iter + 1 << "\n";
            break;
        }
        if (iter == MAX_ITER - 1) {
            std::cerr << "Аварийный выход: достигнуто максимальное число итераций!\n";
        }
    }

    // Запись результата
    std::ofstream result(output_file);
    for (const auto& xi : x) result << xi << "\n";
}

int main() {
    setlocale(LC_ALL, "Russian");
    double omega;
    std::ifstream config("config.txt");
    if (!config.is_open()) {
        std::cerr << "Ошибка открытия файла конфигурации!\n";
        return 1;
    }
    config >> omega;

    // тумблер
    solve_SLAE("config.txt", "matrix.txt", "rhs.txt", "initial_guess.txt", "result.txt", false, omega);  // Метод Якоби
    solve_SLAE("config.txt", "matrix.txt", "rhs.txt", "initial_guess.txt", "result.txt", true, omega);   // Метод Гаусса-Зейделя
    return 0;
}
