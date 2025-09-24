#include <iostream>
// Array Yang Dinamis
#include <vector>
#include <cmath>
using namespace std;

// Mode 0 berdasarakn Error Treshold
// Mode 1 berdasarakn Iterasi
void gaussSeidel(vector<vector<double>> &persIn, vector<double> hasil, short *jmlpers, short *var, bool mode)
{
    int iterasi = 0;
    double err = 100, errTres = 0.0, temp = 0.0, pembagi = 0.0, tempHasil[*var - 1];
    vector<double> hasilAwal;
    cout << hasil[0] << " " << hasil[1] << " " << hasil[2] << endl;
    if (!mode)
    {
        cout << "Err Treshold: ";
        cin >> errTres;
        while (err > errTres)
        {
            hasilAwal = hasil;
            cout << ">> Iterasi ke - " << iterasi << endl;
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
                hasil[i] = temp / pembagi;
                cout << "Pers ke - (" << i + 1 << ") Err= " << abs(tempHasil[i] - hasil[i]) << " , hasil : " << hasil[i] << endl;
                cout << temp << " << Temp || pembagi >> " << pembagi << endl;
                if (i == *jmlpers - 1)
                {
                    for (int h = 0; h < *jmlpers; h++)
                    {
                        err += abs(tempHasil[h] - hasil[h]);
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

void gaussJacobi(vector<vector<double>> &persIn, vector<double> hasil, short *jmlpers, short *var, bool mode)
{
    int iterasi = 0;
    double err = 100, errTres = 0.0, temp = 0.0, pembagi = 0.0, tempHasil[*var - 1];
    cout << hasil[0] << " " << hasil[1] << " " << hasil[2] << endl;
    if (!mode)
    {
        cout << "Err Treshold: ";
        cin >> errTres;
        while (err > errTres)
        {
            cout << ">> Iterasi ke - " << iterasi << endl;
            for (int i = 0; i < *jmlpers; i++)
            {
                temp = 0;
                tempHasil[i] = hasil[i];
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
                    temp += -(persIn[i][j]) * hasil[j];
                }
                hasil[i] = temp / pembagi;
                cout << "Pers ke - (" << i + 1 << ") Err= " << abs(tempHasil[i] - hasil[i]) << " , hasil : " << hasil[i] << endl;
                cout << temp << " << Temp || pembagi >> " << pembagi << endl;
                if (i == *jmlpers - 1)
                {
                    for (int h = 0; h < *jmlpers; h++)
                    {
                        err += abs(tempHasil[h] - hasil[h]);
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

// Pembuat Warna di Terminal
string color(const string &text, const string &code)
{
    return "\033[" + code + "m" + text + "\033[0m";
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
        // Untuk Nambahin Var buat hasil
        var = var + 1;

        if (cin.fail())
        {
            string err = "Masukan Tidak Valid!", code = "1;31"; // merah
            cout << endl
                 << color(err, code) << endl;

            cin.clear();
            cin.ignore(1000, '\n');
        }
        else
        {
            break;
        }
    }

    // Menanyakan Mau berapa Persamaan
    while (true)
    {
        cout << "\n>> Mau Pakai Berapa Pers?: ";
        cin >> jmlPers;

        if (cin.fail())
        {
            string err = "Masukan Tidak Valid!", code = "1;31"; // merah
            cout << endl
                 << color(err, code) << endl;

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
                    string err = "Rumus Tidak Valid!", code = "1;31"; // merah
                    cout << endl
                         << color(err, code) << endl;

                    cin.clear();
                    cin.ignore(1000, '\n');
                }
                else
                {
                    break;
                }
            }
        }
    }
    cout << ">> Kondisi Awal Var?: " << endl;
    for (int i = 0; i < var - 1; i++)
    {
        while (true)
        {
            cout << "X" << i + 1 << " : ";
            cin >> hasil[i];

            if (cin.fail())
            {
                string err = "Masukan Tidak Valid!", code = "1;31"; // merah
                cout << endl
                     << color(err, code) << endl;

                cin.clear();
                cin.ignore(1000, '\n');
            }
            else
            {
                break;
            }
        }
    }
    // Menjadi pengali 1
    hasil[var - 1] = 1;
    // Loop Menu
    while (true)
    {
        cout << "Persamaan yang Mau dihitung ada [" << jmlPers << "] Dengan masing-masing memiliki [" << var << "] variabel" << endl;
        for (int i = 0; i < jmlPers; i++)
        {
            cout << "Persamaan (" << i + 1 << ") ";
            for (int j = 0; j < var; j++)
            {
                // Menampilkan semua persamaan kecuali nilai Z
                if (j < var - 1)
                {
                    if (j > 0)
                    {
                        if (persIn[i][j] >= 0)
                        {
                            cout << "+" << persIn[i][j] << "X" << j + 1 << " ";
                            continue;
                        }
                        // karena kalau dia minus langsung ada minusnya di terminal
                    }
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
                string err = "Masukan Tidak Valid!", code = "1;31"; // merah
                cout << endl
                     << color(err, code) << endl;
                cin.clear();
                cin.ignore(1000, '\n');
            }
            else if (menu > 3)
            {
                string err = "Out Of Index!", code = "1;31"; // merah
                cout << endl
                     << color(err, code) << endl;
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
    vector<vector<double>> pers;
    vector<double> hasil;

    mainMenu(pers, hasil);
    return 0;
}
