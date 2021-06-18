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

```markdown
for (double n = N; n >= 0; n--) {
            A += (1/factorial(2*n+1));
            if (n == 0)
                A *= x;
            else
                A *= -(x*x);
```
Hasil Output :
( Dengan N = 100 )
A : -0.9879136391243737
Time : 0.001511368 detik

( Dengan N = 500 )
A : -0.9879136391243737
Time : 0.005108622 detik

( Dengan N = 1000 )
A : -0.9879136391243737
Time : 0.008196645 detik

( Dengan N = 5000 )
A : -0.9879136391243737
Time : 0.106724261 detik

( Dengan N = 10000 )
A : -0.9879136391243737
Time : 0.415457556 detik

 ( Dengan N = 50000 )
A : -0.9879136391243737
Time : 9.187904015 detik

Dari hasil di atas ternyata membawa hasil yang berbeda dari kedua teknik di atas. Hal ini dapat dilihat bahwa pada saat nilai N = 500, pada model program (a) menampilkan nilai NaN alias tak definisi. Tetapi, untuk model program (b) akan muncul nilainya. Hal ini terjadi karena proses perhitungan pada model program (a) kemungkinan ada nilai yang infinity sehingga hasil akhirnya adalah NaN pada saat perhitungan nilai dengan Math.pow dan method factorial. Untuk model program (b) karena menggunakan teknik pengurangan operasi, maka hasil outputnya bisa tampilkan bahkan pada saat N = 50000. Jika dilihat, model program (a) sangat tidak efisien daripada model program (b) karena pada model (a) operasinya sangat banyak dan membutuhkan waktu yang banyak untuk masuk method Math dan factorial.
Kita bisa lihat bahwa dalam pemrograman, perhitungan data yang nilai sangat besar seperti besar pertumbuhan ekonomi di Indonesia ini, besar energi radiasi matahari atau nilai yang sangat kecil seperti jari – jari atom, besar muatan listrik pada elektron, dan lain – lain sangat berpengaruh dalam program komputer. Tentunya, tipe data yang digunakan harus tepat. Misalkan jika perhitungan data aritmatika biasa seperti deret aritmatika cukup pakai INTEGER, menghitung nilai rata – rata mahasiswa menggunakan FLOAT karena bisa menghasilkan nilai desimal. Oleh karena itu, perlu diperhatikan tipe data-nya beserta bentuk formula yang dibuat seperti menghitung nilai cosinus dengan teknik pengurangan operasi.
                Bagaimana pun juga, metode numerik ini sangat berguna untuk penerapan dalam kehidupan sehari – hari terutama dalam dunia teknik dan sains. Contoh yang sangat sederhana adalah ketika ingin memotong kayu dengan panjang 10 cm, saat mengukur ukuran harus diberi error sekitar 0,1 cm. Hal ini bertujuan agar saat memotong kayu akan menghasilkan ukuran pas 10 cm dan 0,1 cm ini akan hilang akibat gesekan / panas yang menyebabkan menjadi aus. Jadi, bisa dikatakan bahwa numerik ini berpengaruh dalam berbagai macam faktor. Misalkan pembangunan gedung tinggi harus memperhatikan pondasi-nya dan tingkat keretakan bahan yang dipakai ( lebih ke arah elastisitas bahan bangunan ).
                Okay sekian dulu ya blog aku di sini. Untuk metode ini masih sangat banyak. Ya saya usahakan saya posting materi selanjutnya. Udah ya sekian dulu aja. Jangan bosan dengan blog aku siap tau bisa berguna buat kalian wkwk... Okay udah dulu ya. Jangan lupa like blog ku ya. Mantapp (y).
                
### Referensi

[https://docplayer.info/32828129-Galat-dalam-komputasi-numerik.html](https://docplayer.info/32828129-Galat-dalam-komputasi-numerik.html)
[http://nurun.lecturer.uin-malang.ac.id/wp-content/uploads/sites/7/2011/09/Materi-9-Struktur-Atom-lanjutan.pdf](http://nurun.lecturer.uin-malang.ac.id/wp-content/uploads/sites/7/2011/09/Materi-9-Struktur-Atom-lanjutan.pdf)
[https://campus-siskom.blogspot.com/2018/02/metode-numerik-error-asal-dan.html](https://campus-siskom.blogspot.com/2018/02/metode-numerik-error-asal-dan.html)
 

# ~ Numerical Solutions of Algebraic and Transcendental Equestion

 Merupakan persamaan yang berisi fungsi transendental dari variabel yang diselesaikan. Persamaan seperti itu sering kali tidak memiliki solusi bentuk tertutup. Contohnya :
 
 <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/3d65e56f9ce8a46059d20b8244c0adbbb4aa051c" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -4.005ex; width:11.013ex; height:9.176ex;" alt="{\displaystyle {\begin{aligned}x&amp;=e^{-x}\\x&amp;=\cos x\\2^{x}&amp;=x^{2}\end{aligned}}}">
 
 Di dalam kerja ilmiah dan teknik sering dijumpai suatu masalah berkenaan dengan upaya menyelesaikan persamaan yang berbentuk :
 
 <img class="alignnone size-full wp-image-71" alt="rumus a1" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-a1.png" width="58" height="29">
 <a href="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-a1.png"><img class="alignnone size-full wp-image-71" alt="rumus a1" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-a1.png" width="58" height="29"></a>
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; persamaan (1.1)
 
 Menyelesaikan persamaan (1.1) maksudnya adalah mencari suatu nilai berkenaan dengan peubah x sedemikian hingga persamaan tersebut bernilai benar. Nilai-nilai yang dimaksud biasanya disebut dengan nilai-nilai akar.

Bila f(x) berbentuk fungsi polinom sederhana (kuadrat, pangkat tiga, atau pangkat empat) maka ada rumus-rumus aljabar (misalnya metode faktorisasi dan metode pembagian suku banyak), ang dapat digunakan untuk menentukan nilai-nilai akarnya. Sebaliknya, bila suatu polinom berderajat lebih tinggi atau berbentuk transenden seperti,

<img style="-webkit-user-select: none;margin: auto;cursor: zoom-in;background-color: hsl(0, 0%, 90%);transition: background-color 300ms;" src="http://blog.ub.ac.id/aldipradana/files/2013/09/a1.jpg" width="811" height="62">

dan seterusnya, tidak tersedia metode aljabar untuk solusinya. Oleh karena itu harus ditempuh dengan cara aproksimasi.
Dalam bagian ini, akan dibicarakan beberapa metode numerik untuk menyelesaikan permasalahan (1.1) dengan f(x)  adalah fungsi aljabar dan/atau transenden.

### METODE BISEKSI (BISECTION METHOD)
 
Dinamakan metode biseksi (Bi Section) didasarkan atas teknis metode ini adalah “belah dua”. Metode Biseksi diformulasikan berdasarkan Teorema konsep dasar kalkulus : nilai antara dan deret Taylor Teorema 1.1 yang menyatakan bahwa : “Bila fungsi f(x) kontinu dalam selang/interval (a,b) dengan f(a) dan f(b) berlawanan tanda, maka f(α)=0 untuk suatu bilangan α sedemikian hingga a<α<b.

Dinamakan metode biseksi (Bi Section) didasarkan atas teknis metode ini adalah “belah dua”. Metode Biseksi diformulasikan berdasarkan Teorema konsep dasar kalkulus : nilai antara dan deret Taylor Teorema 1.1 yang menyatakan bahwa : “Bila fungsi f(x) kontinu dalam selang/interval (a,b) dengan f(a) dan f(b) berlawanan tanda, maka f(α)=0 untuk suatu bilangan α sedemikian hingga a<α<b.

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/3.png" width="69" height="42">

Bila <img class="alignnone size-full wp-image-77" alt="rumus 4b" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b1.png" width="53" height="18"> atau <img class="alignnone size-full wp-image-78" alt="rumus 4b2" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b2.png" width="32" height="17"> mendekati nilai 0 untuk suatu nilai toleransi yang diberikan maka x0 adalah nilai akar dari <img class="alignnone size-full wp-image-77" alt="rumus 4b" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b1.png" width="53" height="18">

Sebaliknya bila <img class="alignnone size-full wp-image-77" alt="rumus 4b" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b1.png" width="53" height="18"> tidak sama dengan 0 atau <img class="alignnone size-full wp-image-78" alt="rumus 4b2" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b2.png" width="32" height="17"> mendekati nilai 0 tetapi **tidak memenuhi suatu nilai toleransi yang diberikan**, maka berdasarkan Teorema 1.1 ada dua kemungkinan yakni nilai akar berada di antara a dan x0 atau nilai akar berada di antara x0 dan b. Dari salah satu kemungkinan ini, metode Biseksi kembali akan digunakan. Secara geometris, metode Biseksi yang dikemukan di atas diilustrasikan melalui gambar grafik berikut ini.

<img class="alignnone size-full wp-image-79" alt="rumus 5" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-5.png" width="396" height="174">

### METODE NEWTON

Metode NEWTON didasarkan pada aproksimasi linear fungsi dan menggunakan prinsip kemiringan (Tangen) kurvanya.

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/12.png" width="497" height="241">

Kalkulasi dengan metode Newton diawali dengan yang tidak terlalu jauh dari sebuah akar, bergerak sepanjang garis linear (kemiringan atau tangen garis) ke perpotongannya di sumbu-x, dan mengambilnya sebagai titik aproksimasi untuk yang berikutnya. Perlakuan ini diteruskan hingga nilai-nilai x dirasakan sukses cukup dekat ke fungsi bernilai nol. Skema kalkulasinya mengikuti segitiga yang dibangun dengan sudut inklinasi dari kemiringan garis pada kurva di yaitu

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/13.png" width="393" height="57">

Aproksimasi berikutnya diteruskan dengan menghitung  x2 dengan skema yang sama dimana nilai  x3 digantikan  oleh x1. Secara umum metode Newton dirumuskan oleh skema berikut ini:

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/14.png" width="235" height="54">

### METODE POSISI SALAH (REGULA FALSI)

Untuk menghitung nilai akar dari f(x)=0 dapat digunakan metode Posisi salah/Regulasi Falsi. Aturan dari metode ini diterangkan secara geometri dalam Gambar. Dalam gambar yang dimaksud, sketsa grafik dari kurva dinyatakan oleh persamaan y=f(x). Akar dari f(x)=0 yang dicari dinyatakan oleh koordinat x dari titik P yang merupakan perpotongan dari kurva y=f(x) dengan sumbu x. Untuk menggunakan aturan RF, diperlukan dua titik A(α,f (α)) dan B0 (x0,f(x0)).

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/211-300x148.png" width="467" height="174">

Proses selanjutnya adalah menghitung nilai x melalui persamaan garis AB0 yang memotong sumbu x di titik P1. Setelah itu dengan menggunakan koordinat titik P1 yakni (x1, 0) dapat ditentukan titik B1 dengan koordinat (x1, f(x1)) . Dengan demikian garis AB1 akan memotong sumbu x di titik P2 dengan koordinat (x2, 0). Demikian proses ini terus dilakukan hingga diperoleh kondisi Pn sangat dekat ke P yakni |Pn-P|<toleransi. Dari proses pencapaian nilai akar di titik P, dihasilkan barisan nilai-nilai x0,x1,x2,……,xn yang diharapkan akan konvergen ke absis x pada titik P, yaitu akar yang dicari dari persamaan yang diberikan.

Persamaan garis AB0 adalah :

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/22.png" width="244" height="61">

Karena persamaan garis AB0 melalui titik P1(x1,0), maka diperoleh :

<img class="alignnone size-medium wp-image-82" alt="a" src="http://blog.ub.ac.id/aldipradana/files/2013/09/a-300x145.jpg" width="243" height="116">

Sumber : [http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&cad=rja&ved=0CEUQFjAD&url=http%3A%2F%2Fmatematikaindo.
files.wordpress.com%2F2010%2F04%2Fmetode-numerik-buku-ajar-unila.pdf&ei=cIA7UrhFi9KtB4-CgcgK&usg=AFQjCNGex_IcMyn18Q6O-ZmCHZ_E3638CQ&sig2=1C2ZFEGmJtwWu
zqhyEtFnQ](http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&cad=rja&ved=0CEUQFjAD&url=http%3A%2F%2Fmatematikaindo.files.wordpress.com%2F2010%2F04%2Fmetode-numerik-buku-ajar-unila.pdf&ei=cIA7UrhFi9KtB4-CgcgK&usg=AFQjCNGex_IcMyn18Q6O-ZmCHZ_E3638CQ&sig2=1C2ZFEGmJtwWuzqhyEtFnQ)








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


