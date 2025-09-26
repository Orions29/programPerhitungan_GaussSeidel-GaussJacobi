#include <iostream>
// Array Yang Dinamis
#include <vector>
#include <iomanip>
using namespace std;

// void gaussSeidelTemplate(vector<vector<double>> &persIn, vector<double> &hasil, short *jmlpers, short *var, double &err, double &errTres, double &pembagi, double *&tempHasilPers)
// {
//     int iterasi = 0;
//     double err = 100, errTres = 0.0, pembagi = 0.0, tempHasilPers[*var - 1];
//     cout << "\033[1;5m>> Iterasi ke - \033[0m" << iterasi << endl;
//     for (int i = 0; i < *jmlpers; i++)
//     {
//         hasil[*var - 1] = 0;
//         tempHasilPers[i] = hasil[i];
//         for (int j = 0; j < *var; j++)
//         {
//             if (j == i)
//             {
//                 pembagi = persIn[i][j];
//                 continue;
//             }
//             else if (j == (*var - 1))
//             {
//                 hasil[j] += persIn[i][j];
//                 continue;
//             }
//             else
//                 hasil[*var - 1] += -(persIn[i][j]) * hasil[j];
//         }
//         hasil[i] = hasil[*var - 1] / pembagi;
//         cout << "Pers ke - (" << i + 1 << ") Hasil = " << hasil[i] << ", \033[31mError Absolute: \033[0m" << abs(tempHasilPers[i] - hasil[i]) << endl;
//         // Buat Debug
//         // cout << temp << " << Temp || pembagi >> " << pembagi << endl;
//         if (i == *jmlpers - 1)
//         {
//             for (int h = 0; h < *jmlpers; h++)
//             {
//                 err += abs(tempHasilPers[h] - hasil[h]);
//             }
//             // Menjadi Error Threshold Rata rata
//             err /= *var - 1;
//             cout << "Rata-rata Error Absolute : " << err << endl;
//         }
//     }
// }

/**
  @brief Fungsi Untuk Penyelesaian non-analitik dengan Metode Gauss Jacobi
  @param persIn Vector Persamaan yang mau diselesaikan
  @param hasil Vector untuk menyimpan hasil perhitungan
  @param jmlpers jumlah persamaan yang mau diselesaikan
  @param var Jumlah variabel yang ada pada masing masing
                pertamaan (Variabel terakhir sebagai konstanta pengali Z_result)
  @param mode - `0` or`False`berhenti di error_value tertentu
                - `1` or `True` berhenti di iterasi
 */
void gaussSeidel(vector<vector<double>> &persIn, vector<double> hasil, short *jmlpers, short *var, bool mode)
{
    // Membersihkan Seluruh Program Sebelumnya
    system("cls");
    int iterasi = 0;
    double err = 100, errTres = 0.0, pembagi = 0.0, tempHasilPers[*var - 1];
    cout << "\033[34m<<Metode Gauss-Seidel>>\033[0m" << endl;
    cout << "\033[32mNilai Awalan:\033[0m \n X1: " << hasil[0] << ", X2: " << hasil[1] << ", X3: " << hasil[2] << endl;
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

        cout << "Operasi Penyelesaian: " << setprecision(10) << endl;
        while (err > errTres)
        {
            cout << "\033[1;5m>> Iterasi ke - \033[0m" << iterasi << endl;
            for (int i = 0; i < *jmlpers; i++)
            {
                hasil[*var - 1] = 0;
                tempHasilPers[i] = hasil[i];
                for (int j = 0; j < *var; j++)
                {
                    if (j == i)
                    {
                        pembagi = persIn[i][j];
                        continue;
                    }
                    else if (j == (*var - 1))
                    {
                        hasil[j] += persIn[i][j];
                        continue;
                    }
                    else
                        hasil[*var - 1] += -(persIn[i][j]) * hasil[j];
                }
                hasil[i] = hasil[*var - 1] / pembagi;
                cout << "Pers ke - (" << i + 1 << ") Hasil = " << hasil[i] << ", \033[31mError Absolute: \033[0m" << abs(tempHasilPers[i] - hasil[i]) << endl;
                // Buat Debug
                // cout << temp << " << Temp || pembagi >> " << pembagi << endl;
                if (i == *jmlpers - 1)
                {
                    for (int h = 0; h < *jmlpers; h++)
                    {
                        err += abs(tempHasilPers[h] - hasil[h]);
                    }
                    // Menjadi Error Threshold Rata rata
                    err /= *var - 1;
                    cout << "Rata-rata Error Absolute : " << err << endl;
                }
            }
            iterasi++;
            // persIn[0]
        }
    }
}

// FIXME - Bentar ini variabel nya baru kuubah ubah belum jadi ini
void gaussJacobi(vector<vector<double>> &persIn, vector<double> hasilPers, short *jmlpers, short *var, bool mode)
{
    int iterasi = 0;
    double err = 100, errTres = 0, temp = 0, pembagi = 0, tempHasil[*var - 1] = {};
    vector<double> hasilAwal;
    cout << hasilPers[0] << " " << hasilPers[1] << " " << hasilPers[2] << endl;
    if (!mode)
    {
        cout << "\033[1;31m Masukkan Err Treshold>>\033[0m ";
        cin >> errTres;
        while (err > errTres)
        {
            hasilAwal = hasilPers;
            cout << ">> Iterasi ke - " << iterasi << "\n";
            for (int i = 0; i < *jmlpers; i++)
            {
                temp = 0;
                tempHasil[i] = hasilAwal[i];
                for (int j = 0; j < *var; j++)
                {
                    if (j == i)
                    {
                        pembagi = persIn[i][j];
                        continue;
                    }
                    else if (j == *var - 1)
                    {
                        temp += persIn[i][j];
                        continue;
                    }
                    temp += -(persIn[i][j]) * hasilAwal[j];
                }
                // Hasil Persamaan
                hasilPers[i] = temp / pembagi;
                cout << "Pers ke - (" << i + 1 << ") Err= " << abs(tempHasil[i] - hasilPers[i]) << " , hasil : " << hasilPers[i] << endl;
                cout << temp << " << Temp || pembagi >> " << pembagi << endl;
                if (i == *jmlpers - 1)
                {
                    for (int h = 0; h < *jmlpers; h++)
                    {
                        err += abs(tempHasil[h] - hasilPers[h]);
                    }
                    // Menjadi Error Threshold Rata rata
                    err /= *var - 1;
                }
            }
            iterasi++;
            // persIn[0]
        }
    }
}

void mainMenu(vector<vector<double>> &persIn, vector<double> &hasil)
{
    short var, jmlPers, menu;
    cout << "\033[44;1;37m<==== Menu Utama Pendekatan Numerik ====>\033[0m" << "\n  => Rio Meidi A / 123240175" << endl
         << endl;
    // Menanyakan mau Berapa Variabel
    while (true)
    {
        cout << "\n>> Mau Pakai Berapa Variabel?: ";
        cin >> var;
        if (cin.fail())
        {
            cout << "\033[1;32mMasukan Tidak Valid!\033[0m";

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
            cout << "\033[1;32mMasukan Tidak Valid!\033[0m";

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
                    cout << "\033[1;32mKonstanta Tidak Valid!\033[0m";
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
                cout << "\033[1;32mMasukan Tidak Valid!\033[0m";

                cin.clear();
                cin.ignore(1000, '\n');
            }
            else
            {
                break;
            }
        }
    }
    // Loop Menu
    while (true)
    {
        cout << "Persamaan yang Mau dihitung ada [" << jmlPers << "] Dengan masing-masing memiliki [" << var - 1 << "] variabel" << endl;
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
        cout << "\t1.) Gauss-Jacobi\n\t2.) Gauss-Seidel\n\t3.) Exit" << endl;
        while (true)
        {
            cout << "\n>> Menu Keberapa: ";
            cin >> menu;

            if (cin.fail())
            {
                cout << "\033[1;32mMasukan Tidak Valid!\033[0m";
                cin.clear();
                cin.ignore(1000, '\n');
            }
            else if (menu > 3)
            {
                cout << "\033[1;32mOut Of Index!\033[0m";
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
            gaussJacobi(persIn, hasil, &jmlPers, &var, 0);
            break;
        case 2:
            gaussSeidel(persIn, hasil, &jmlPers, &var, 0);
            break;
        case 3:
            exit(0);
            break;
        default:
            break;
        }
    }
}

int main(int argc, char const *argv[])
{
    // Inisialisasi Tempat Untuk Persamaan dan Hasil
    vector<vector<double>> pers = {};
    vector<double> hasil = {};

    mainMenu(pers, hasil);
    return 0;
}
