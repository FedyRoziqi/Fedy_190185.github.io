# Komputasi Numerik
# ~ Errors in Numerical Computations
Atau yang biasa kita sebut sebagai nilai eror,lalu apasih nilai eror itu :

Secara umum,terdapat dua sumber utama penyebab terjadinya error dalam perhitungan numerik, yaitu:

### 1. Error pembulatan ( round-off error )
Perhitungan dengan metode numerik hampir selalu menggunakan bilangan riil.Masalah timbul apabila komputasi numerik dikerjakan oleh mesin (dalam hal ini dengan menggunakan komputer) karena semua bilangan riil tidak dapat disajikan secara tepat di dalam komputer. Keterbatasan komputer dalam menyajikan bilangan riil menghasilkan error yang disebut error pembulatan.

Sebagai contoh 1/6 = 0.166666666… tidak dapat dinyatakan secara tepat oleh komputer karena digit 6 panjangnya tidak terbatas. Komputer hanya mampu merepresentasikan sejumlah digit (atau bit dalam sistem biner) saja. Bilangan riil yang panjangnya melebihi jumlah digit (bit) yang dapat direpresentasikan oleh komputer dibulatkan ke bilangan terdekat. Misalnya sebuah komputer hanya dapat merepresentasikan bilangan riil dalam 6 digit angka berarti, maka representasi bilangan 1/6 = 0.1666666666… di dalam komputer 6-digit tersebut adalah 0.166667. Galat pembulatannya adalah 1/6 – 0.166667 = -0.000000333.

Contoh dalam sistem biner misalnya 1/10 = 0.000110011001100110011 00110011…2 direpresentasikan di dalam komputer dalam jumlah bit yang terbatas. Kebanyakan komputer digital mempunyai dua buah cara penyajian bilangan riil, yaitu bilangan titik-tetap (fixed point) dan bilangan titik-kambang (floatingpoint). Dalam format bilangan titik -tetap setiap bilangan disajikan dengan jumlah tempat desimal yang tetap, misalnya 62.358, 0.013, 1.000. Sedangkan dalam format bilangan titik-kambang setiap bilangan disajikan dengan jumlah digit berarti yang sudah tetap, misalnya 0.6238 103 0.1714 ^10-13 atau ditulis juga 0.6238E+03 0.1714E-13. Digit-digit berarti di dalam format bilangan titik-kambang disebut juga angka bena (significant figure).

### 2. Galat Pemotongan ( truncation error )

Galat pemotongan adalah galat yang ditimbulkan oleh pembatasan jumlah komputasi yang digunakan pada proses metode numerik. Banyak metode dalam metode numerik yang penurunan rumusnya menggunakan proses iterasi yang jumlahnya tak terhingga, sehingga untuk membatasi proses penghitungan, jumlah iterasi dibatasi sampai langkah ke n. Hasil penghitungan sampai langkah ke n akan menjadi hasil hampiran dan nilai penghitungan langkah n keatas akan menjadi galat pemotongan. dalam hal ini galat pemotongan kan menjadi sangat kecil sekali jika nilai n di perbesar. Konsekuensinya tentu saja jumlah proses penghitungannya akan semakin banyak.

### 3. Cara atau penyebab lainnya antara lain :
#### 1. Round off errors
  Kesalahan ini biasanya akibat proses pembulatan dalam perhitungan. Secara umum, proses pembulatan ada 2 aturan yaitu :
-          Jika digit yang dibulatkan kurang dari 5, maka tidak terjadi pembulatan.
-          Sebaliknya, jika lebih dari 5, maka terjadi pembulatan yaitu dengan menambah satu.

#### 2. Truncation errors
  Kesalahan pemotongan biasanya terjadi karena pembuangan suku yang berderajat tinggi. Sebagai contoh untuk menghitung nilai cosinus dapat menggunakan deret Taylor
  
#### 3.  Range errors
Untuk kesalahan ini berkaitan dengan batas dalam jangkauan representasi angka. Ini bisa dikatakan bahwa jika hasil perhitungan melebihi jangkauan, maka komputer akan menampilkan hasil yang tidak beraturan ( anggap saja hasil yang diperoleh diatur lagi oleh OS yang kita pakai ).
-          Menghitung jari – jari atom Bohr.
-          Menghitung nilai sinus dengan teknik Taylor.
Bentuk rumus-nya sudah ditulis sebelumnya. Disini akan diberi 2 teknik perhitungan dalam penggunaan deret Taylor.

a.       Teknik Penggunaan Fungsi Math.
Bentuk program dapat dibuat seperti di bawah ini :

```markdown
for (double n = 0; n <= N; n++)
            B += (Math.pow(-1, n))*(Math.pow(x, 2 * n+1))/(factorial(2*n+1));
```
 Hasil output :
( Dengan N = 100 )
A : -0.9881209821429402
Time : 7.33534E-4 detik

( Dengan N = 500 )
A : NaN
Time : 0.002157699 detik

( Dengan N = 1000 )
A : NaN
Time : 0.005704589 detik

( Dengan N = 5000 )
A : NaN
Time : 0.108811078 detik

( Dengan N = 10000 )
A : NaN
Time : 0.418405214 detik

( Dengan N = 50000 )
A : NaN
Time : 24.734975616 detik

b.      Teknik Pengurangan Operasi.


 




You can use the [editor on GitHub](https://github.com/FedyRoziqi/fedyrz.github.io/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/FedyRoziqi/fedyrz.github.io/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.


