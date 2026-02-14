#include <iostream>
#include <cmath>
#include <iomanip> 
using namespace std;

double ParabolaMethod(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        sum += 2 * f(a + i * h);
    }
    for (int i = 0; i < n; i++) {
        sum += 4 * f(a + i * h + h / 2.0);
    }
    return (h / 6.0) * sum;
}

double TestFunc(double x) {
    return 6 * pow(x, 5);
}

double VariantFunc(double x) {
    return pow(x, 1.0 / 13.0) * sqrt(1 + x * x);
}

void RunTestTask() {
    cout << "\nОтладочный пример функции f(x)=6x^5 на диапазоне [0, 1] \n";
    cout << setw(6) << "n" << " | "
        << setw(10) << "K дельта" << " | "
        << setw(15) << "дельта точное" << " | "
        << setw(15) << "дельта Рунге" << " | "
        << setw(15) << "дельта теор" << endl;
    cout << string(75, '-') << endl;

    double a = 0.0, b = 1.0;
    double exact_value = 1.0; //точное значение

    double prev_integral = 0.0; //пред. значение интеграла
    double prev_error_runge = 0.0; // пред. оценка погрешности

    double M4 = 720.0; 

    for (int n = 1; n <= 65536; n *= 2) {
        double integral = ParabolaMethod(TestFunc, a, b, n);

        // дельта точное
        double delta_exact = exact_value - integral;
        // дельта рунге
        double delta_runge = 0.0;
        if (n > 1) {
            delta_runge = (integral - prev_integral) / (pow(2,4)-1);
        }

        // дельта теор
        double h = (b - a) / n;
        double delta_teor = (M4 / 2880.0) * (b - a) * pow(h, 4);

        // K дельта
        double k_delta = 0.0;
        if (n > 2 && abs(delta_runge) > 1e-20) {
            k_delta = abs(prev_error_runge / delta_runge);
        }

        // Вывод
        cout << setw(6) << n << " | ";
        if (n > 2) cout << setw(10) << fixed << setprecision(1) << k_delta << " | ";
        else       cout << setw(10) << "-" << " | ";

        cout << scientific << setprecision(5)
            << setw(15) << delta_exact << " | ";

        if (n > 1) cout << setw(15) << delta_runge << " | ";
        else       cout << setw(15) << "-" << " | ";

        cout << setw(15) << delta_teor << endl;

        prev_integral = integral;
        prev_error_runge = delta_runge;
    }
}

void RunVariantTask(double a, double b) {
    cout << "\nИнтегрирование функции x^(1/13)*(1+x^2)^(1/2) на отрезке [" << a << ", " << b << "] \n";

    cout << setw(6) << "n" << " | " << setw(14) << "интеграл" << " | " << setw(20) << "дельта Рунге" << " | " << setw(10) << "K дельта" << endl;
    cout << string(80, '-') << endl;

    double epsilon = 1e-6;
    int n = 1;

    double prev_integral = ParabolaMethod(VariantFunc, a, b, n);
    double integral = 0.0;
    double prev_error = 0.0;
    cout << setw(6) << n << " | " << fixed << setprecision(12) << prev_integral << " | "
        << scientific << setprecision(5) << setw(20) << "-" << " | " << setw(10) << "-" << endl;

    while (true) {
        n *= 2;
        integral = ParabolaMethod(VariantFunc, a, b, n);

        double error_runge = (integral - prev_integral) / (pow(2, 4) - 1);

        double k_delta = 0.0;
        
        if (abs(prev_error) > 1e-20) {
            k_delta = abs(prev_error / error_runge);
        }

        cout << defaultfloat;
        cout << setw(6) << n << " | " << fixed << setprecision(12) << integral << " | "
            << scientific << setprecision(5) << setw(20) << error_runge << " | ";

        if (abs(prev_error) > 1e-20)
            cout << fixed << setprecision(1) << setw(10) << k_delta << endl;
        else
            cout << setw(10) << "-" << endl;

        if (abs(error_runge) < epsilon) {
            cout << string(65, '-') << endl;
            cout << "Точность " << scientific << epsilon << " достигнута при n = " << n << endl;
            cout << "Результат: " << fixed << setprecision(9) << integral << endl;
            break;
        }

        prev_integral = integral;
        prev_error = error_runge;

        if (n > 2000000) break;
    }
}
int main() {
    setlocale(LC_ALL, "ru");
    RunTestTask();
    RunVariantTask(0.0, 1.5);
    RunVariantTask(0.001, 1.5);
    return 0;
}