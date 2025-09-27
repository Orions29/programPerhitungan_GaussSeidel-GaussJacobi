/**
 * @file persamaan.cpp
 * @author Rio Meidi Ataurrahman (riomeidi29@gmail.com)
 * @brief Untuk Tugas Komputer Dan Numerasi ke 1
 * @version 0.1
 * @date 2025-09-27
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <iostream>
// Array Yang Dinamis
#include <vector>
#include <iomanip>
using namespace std;

/**
 * @brief Template Engine untuk menjalankan perhitungan Gauss-Jacobi dan Gauss-Seidel
 *
 * @param persIn Array Persamaan yang berisi konstanta dan nilai Z pada akhir
 * @param hasil Array untuk menyimpan hasil persamaan
 * @param hasilPrev Array untuk menyimpan hasil sebelum perhitungan
 * @param jmlpers Jumlah Persamaan yang mau dihitung
 * @param var Jumlah Variabel pada tiap Persamaan
 * @param err Mean Absolute Error pada tiap hasil persamaan sebelum dan setelah perhitungan
 * @param iterasi jumlah iterasi yang terlah dilakukan
 * @param mode `0` / `False` untuk Gauss-Jacobi dan \enmd `1` / `True` untuk Gauss-Seidel
 */
void gaussSeidelTemplate(const vector<vector<double>> &persIn, vector<double> &hasil, vector<double> &hasilPrev, const short *jmlpers, const short *var, double &err, int &iterasi, bool mode)
{
    hasilPrev = hasil;
    double pembagi = 0.0, tempErr = 0.0;
    err = 0;
    cout << "\033[1;5m>> Iterasi ke - \033[0m" << iterasi << endl;
    for (int i = 0; i < *jmlpers; i++)
    {
        tempErr = 0;
        hasil[*var - 1] = 0;
        pembagi = persIn[i][i];
        for (int j = 0; j < *var - 1; j++)
        {
            if (j != i)
            {
                if (!mode)
                {
                    hasil[*var - 1] += -(persIn[i][j]) * hasilPrev[j];
                }
                else
                    hasil[*var - 1] += -(persIn[i][j]) * hasil[j];
            }
        }
        // Hasil Persamaan
        hasil[i] = (hasil[*var - 1] + persIn[i][*var - 1]) / pembagi;
        tempErr = abs(hasilPrev[i] - hasil[i]);
        err += tempErr;
        cout << "Pers ke - (" << i + 1 << ") Hasil = " << hasil[i] << ", \033[31mError Absolute: \033[0m" << tempErr << endl;
        // Buat Debug
        // cout << temp << " << Temp || pembagi >> " << pembagi << endl;
    }
    // Menjadi Error Rata rata
    err /= *jmlpers;
    cout << "Rata-rata Error Absolute : " << err << endl;
    iterasi++;
}

/**
 *  @brief Fungsi Untuk Penyelesaian non-analitik dengan Metode Gauss Jacobi
 *  @param persIn Vector Persamaan yang mau diselesaikan
 *  @param hasil Vector untuk menyimpan hasil perhitungan
 *  @param jmlpers jumlah persamaan yang mau diselesaikan
 *  @param var Jumlah variabel yang ada pada masing masing
 *             pertamaan (Variabel terakhir sebagai konstanta pengali Z_result)
 *  @param mode - `0` or`False`berhenti di error_value tertentu
 *              - `1` or `True` berhenti di iterasi
 */
void gaussMenu(const vector<vector<double>> &persIn, vector<double> hasil, const short *jmlpers, const short *var, const bool menu)
{
    // Membersihkan Seluruh Program Sebelumnya
    system("cls");
    bool mode, pil;
    int iterasi = 0, iterasiMax = 0;
    double err = 0, errTres = 0.0, pembagi = 0.0;
    // Buat Gauss jacobi
    vector<double> hasilPrev;
    if (!menu)
    {
        cout << "\033[34m<<Metode Gauss-Jacobi>>\033[0m" << endl;
    }
    else
        cout << "\033[34m<<Metode Gauss-Seidel>>\033[0m" << endl;
    cout << "\033[32mNilai Awalan:\033[0m \n X1: " << hasil[0] << ", X2: " << hasil[1] << ", X3: " << hasil[2] << endl;
    // Error Handling Mode
    do
    {
        cout << "Mode Apa \033[31m(0) Max Err\033[0m \033[32m(1) Max Iterasi :\033[0m ";
        cin >> mode;
        if (cin.fail())
        {
            cout << "\033[1;31mError Input\033[0m - Woops! Inputkan lagi";
        }
        else
            break;
    } while (true);
    if (!mode)
    {
        // Error Handling memasukkan float
        do
        {
            cout << "\033[31mErr Treshold:\033[0m ";
            cin >> errTres;
            if (cin.fail())
            {
                cout << "\033[1;31mError Input\033[0m - Woops! Inputkan lagi";
            }
            else
                break;
        } while (true);

        cout << "Operasi Penyelesaian: " << fixed << setprecision(10) << endl;
        do
        {
            gaussSeidelTemplate(persIn, hasil, hasilPrev, jmlpers, var, err, iterasi, menu);
        } while (err > errTres);
    }
    else
    {
        // Error Handling memasukkan Iterasi ke
        do
        {
            cout << "\033[32mBerapa Perulangan:\033[0m ";
            cin >> iterasiMax;
            if (cin.fail())
            {
                cout << "\033[1;31mError Input\033[0m - Woops! Inputkan lagi";
            }
            else
                break;
        } while (true);
        cout << "Operasi Penyelesaian: " << fixed << setprecision(10) << endl;
        for (int i = 0; i <= iterasiMax; i++)
        {
            gaussSeidelTemplate(persIn, hasil, hasilPrev, jmlpers, var, err, iterasi, menu);
        }
    }
}

/**
 * @brief Mengecek Konvergensi Diagonal dan Mengecek Dominansi Diagonal
 *
 * @param persIn Koefisien persamaan yang mau dicek
 * @param jmlPers Jumlah Persamaan
 * @param var Jumlah Variabel Tiap Persamaan
 */
void checkAndSwapKonvergenDominan(vector<vector<double>> &persIn, const short jmlPers, const short var, string &log)
{
    for (int i = 0; i < jmlPers; ++i)
    {
        double maxKoef = 0.0;
        int maxInRow = i;
        // Akan Memulai Perulangan dimulai dengan baris selanjutnya, Baris Sebelumnya yang sudah dicek akan dibiarkan
        for (int k = i; k < jmlPers; ++k)
        {
            if (abs(persIn[k][i]) > maxKoef)
            {
                maxKoef = abs(persIn[k][i]);
                maxInRow = k;
            }
        }
        if (maxInRow != i)
        {
            swap(persIn[i], persIn[maxInRow]);
            cout << "\033[1;33m Baris " << maxInRow + 1 << " ditukar dengan Baris " << i + 1
                 << " untuk Dominasi Diagonal!\033[0m" << endl;
        }
        // Ngecek apakah persamaan ini dominan atau enggak
        double sumNonDiagonal = 0.0;
        for (int j = 0; j < var - 1; ++j)
        {
            if (i != j)
            {
                sumNonDiagonal += abs(persIn[i][j]);
            }
        }
        if (abs(persIn[i][i]) > sumNonDiagonal)
        {
            log += "\033[1;32m >> Persamaan " + to_string(i + 1) + " Sudah Dominan Diagonal.\033[0m\n";
        }
        else
        {
            log += "\033[1;5;33m >> Persamaan " + to_string(i + 1) + "\033[0m ,Tidak Strictly Dominan Diagonal (Mungkin jadi nih ye tidak konvergen).\n";
        }
    }
}

/**
 * @brief Untuk Menampilkan Menu dan Looping Main Menu
 *
 * @param persIn Persamaan yang mau diselesaikan
 * @param hasil Hasil dari persamaan nantinya
 */
void mainMenu(vector<vector<double>> &persIn, vector<double> &hasil)
{
    bool loopMain = false, loopSubMenu = true;
    do
    {
        system("cls");
        loopSubMenu = true;
        short var, jmlPers, menu;
        string logErr;
        cout << "\033[44;1;37m<==== Menu Utama Pendekatan Numerik ====>\033[0m" << "\n  => Rio Meidi A / 123240175" << endl
             << endl;
        // Menanyakan mau Berapa Variabel
        while (true)
        {
            cout << "\n>> Mau Pakai Berapa Variabel?: ";
            cin >> var;
            if (cin.fail())
            {
                cout << "\033[1;31mMasukan Tidak Valid!\033[0m";

                cin.clear();
                cin.ignore(1000, '\n');
            }
            else
            {
                var++;
                break;
            }
            // Menambahkan Var untuk Konstanta Z
        }
        // Menanyakan Mau berapa Persamaan
        while (true)
        {
            cout << "\n>> Mau Pakai Berapa Pers?: ";
            cin >> jmlPers;

            if (cin.fail())
            {
                cout << "\033[1;31mMasukan Tidak Valid!\033[0m";

                cin.clear();
                cin.ignore(1000, '\n');
            }
            else
            {
                break;
            }
        }
        // Mereserve Memory untuk sejumlah persamaan dan variabelnya
        persIn.resize(jmlPers);
        hasil.resize(var);
        // Resize Tiap Dimensi dari Array Persamaan
        for (int i = 0; i < jmlPers; i++)
        {
            persIn[i].resize(var);
        }
        // Memasukkan Persamaan
        for (int i = 0; i < jmlPers; i++)
        {
            cout << "\n >> Masukkan Persamaan ke - " << i + 1 << endl;
            for (int j = 0; j < var; j++)
            {
                while (true)
                {
                    if (j == (var - 1))
                    {
                        cout << "\n\t> Hasil Z : ";
                        cin >> persIn[i][j];
                    }
                    else
                    {
                        cout << "\n\t> Konstanta X" << j + 1 << " : ";
                        cin >> persIn[i][j];
                    }
                    if (cin.fail())
                    {
                        cout << "\033[1;31mKonstanta Tidak Valid!\033[0m";
                        cin.clear();
                        cin.ignore(1000, '\n');
                        continue;
                    }
                    break;
                }
            }
        }
        cout << ">> Kondisi Awal Var?: " << endl;
        // Var - 1 karena hasil temp persamaan ditaruh di belakang sendiri
        for (int i = 0; i < var - 1; i++)
        {
            while (true)
            {
                cout << "X" << i + 1 << " : ";
                cin >> hasil[i];

                if (cin.fail())
                {
                    cout << "\033[1;31mMasukan Tidak Valid!\033[0m";

                    cin.clear();
                    cin.ignore(1000, '\n');
                }
                else
                {
                    break;
                }
            }
        }
        // Cek apakah dia konvergen atau tidak
        checkAndSwapKonvergenDominan(persIn, jmlPers, var, logErr);
        // Loop Menu
        do
        {
            cout << logErr;
            cout << "Persamaan yang Mau dihitung ada [" << jmlPers << "] Dengan masing-masing memiliki [" << var - 1 << "] variabel" << endl;
            cout.unsetf(ios::fixed);
            for (int i = 0; i < jmlPers; i++)
            {
                cout << "Persamaan (" << i + 1 << ") ";
                for (int j = 0; j < var; j++)
                {
                    // Menampilkan semua persamaan kecuali nilai Z
                    if (j < var - 1)
                    {
                        // Depan sendiri gaperlu dikasi penanda +
                        if (j > 0)
                        {
                            if (persIn[i][j] >= 0)
                            {
                                cout << "+" << persIn[i][j] << "X" << j + 1 << " ";
                                continue;
                            }
                            // karena kalau dia minus langsung ada minusnya di terminal
                        }
                        // Kalau dia depan sendiri tampilkan apa adanya
                        cout << persIn[i][j] << "X" << j + 1 << " ";
                        continue;
                    }
                    cout << "= " << persIn[i][j];
                }
                cout << endl;
            }

            cout << "Pilih Mau Pakai Metode apa?: \n";
            cout << "\t1.) Gauss-Jacobi\n\t2.) Gauss-Seidel\n\t3.) Persamaan Baru\n\t4.) Exit" << endl;
            while (true)
            {
                cout << "\n>> Menu Keberapa: ";
                cin >> menu;

                if (cin.fail())
                {
                    cout << "\033[1;31mMasukan Tidak Valid!\033[0m";
                    cin.clear();
                    cin.ignore(1000, '\n');
                }
                else if (menu > 4)
                {
                    cout << "\033[1;31mOut Of Index!\033[0m";
                    cin.clear();
                    cin.ignore(1000, '\n');
                }
                else
                {
                    break;
                }
            }
            switch (menu)
            {
            case 1:
                gaussMenu(persIn, hasil, &jmlPers, &var, 0);
                break;
            case 2:
                gaussMenu(persIn, hasil, &jmlPers, &var, 1);
                break;
            case 3:
                loopSubMenu = false;
                loopMain = true;
                break;
            case 4:
                loopSubMenu = false;
                loopMain = false;
                break;
            default:
                // barangkali i think..
                cout << "\033[1;5;31mPlis ini Error handling terakhir....\033[0m" << endl;
                system("pause");
                break;
            }
        } while (loopSubMenu);
    } while (loopMain);
}
int main(int argc, char const *argv[])
{
    // Inisialisasi Tempat Untuk Persamaan dan Hasil
    vector<vector<double>> pers = {};
    vector<double> hasil = {};

    mainMenu(pers, hasil);
    return 0;
}
