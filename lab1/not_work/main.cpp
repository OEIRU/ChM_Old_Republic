#include "main.h"

int main(int argc, char** argv)
{
    int memory_num_first, memory_num_end;
    int choice = -1;
    std::string inputFile, outputFile;
    inputFile = "test1.txt";
    while (choice != 7) {
        showMenu();
        std::cout << "Введите номер команды: ";
        std::cin >> choice;

        switch (choice) {
            case 1: {
                std::fstream ifile(inputFile);
                struct SkylineStorageMatrix<double> a;
                readSkyline(ifile, a);
                std::vector<double> b;
                readVector(ifile, b);

                outputFile = "solve_double.txt";
                std::fstream ofile(outputFile, std::ios::out);
                if (!ofile.is_open()) 
                    std::cerr << "Ошибка: не удалось открыть файл " << outputFile << " для записи." << std::endl;
                std::vector<double> x(b.size());
                solve(x, a, b);
                writeVector(x, ofile);
                break;
            }
            case 2: {
                std::fstream ifile(inputFile, std::ios::in);
                struct SkylineStorageMatrix<float> a;
                readSkyline(ifile, a);
                std::vector<float> b;
                readVector(ifile, b);

                outputFile = "solve_float.txt";
                std::fstream ofile(outputFile, std::ios::out);
                if (!ofile.is_open()) 
                    std::cerr << "Ошибка: не удалось открыть файл " << outputFile << " для записи." << std::endl;
                std::vector<float> x(b.size());
                solve(x, a, b);
                writeVector(x, ofile);
                break;
            }
            case 3: {
                std::fstream ifile(inputFile, std::ios::in);
                struct SkylineStorageMatrix<double> a;
                readSkyline(ifile, a);
                std::vector<double> b;
                readVector(ifile, b);
                std::cout << "Отладка: Исходный файл прочитан. " << std::endl;
                std::vector<std::vector<double>> A;
                skylineToSquare(a, A);
                std::cout << "Отладка: Перевод в квадрат переведен. " << std::endl;

                outputFile = "solve_square.txt";    
                std::fstream ofile(outputFile, std::ios::out);
                if (!ofile.is_open()) 
                    std::cerr << "Ошибка: не удалось открыть файл " << outputFile << " для записи." << std::endl;   
                std::vector<double> x(b.size());
                solveGauss(A, x, b);
                std::cout << "Отладка: Гаусс посчитан. " << std::endl;

                writeVector(x, ofile);
                std::cout << "Отладка: Запись результата. " << std::endl;

                break;
            }
            case 4: {
                std::cout << "Введите начальное значение: " << std::endl;
                std::cin >> memory_num_first; 
                std::cout << "Введите конечное значение: " << std::endl;
                std::cin >> memory_num_end;
                hilbertTestSeries<int>(memory_num_first, memory_num_end);
                break;
            }
            case 5: {
                std::cout << "Введите размерность матрицы: " << std::endl;
                std::cin >> memory_num_first; 
                conditionNumberTestSeries(memory_num_first);
                break;
            }
            case 6: {
                std::cout << "Введите размерность матрицы: " << std::endl;
                std::cin >> memory_num_first; 
                compareGaussLusq(memory_num_first);
                break;
            }
            case 7: {
                std::cout << "Выход из программы." << std::endl;
                break;
                
            }
            default: {
                std::cout << "Неверный выбор. Пожалуйста, попробуйте снова." << std::endl;
                break;
            }
        }
    }

    return 0;
}