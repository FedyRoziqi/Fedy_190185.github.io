# Komputasi Numerik
# ~ Errors in Numerical Computations
Atau yang biasa kita sebut sebagai nilai eror,lalu apasih nilai eror itu :

Secara umum,terdapat dua sumber utama penyebab terjadinya error dalam perhitungan numerik, yaitu:

**Error pembulatan ( round-off error )**
Perhitungan dengan metode numerik hampir selalu menggunakan bilangan riil.Masalah timbul apabila komputasi numerik dikerjakan oleh mesin (dalam hal ini dengan menggunakan komputer) karena semua bilangan riil tidak dapat disajikan secara tepat di dalam komputer. Keterbatasan komputer dalam menyajikan bilangan riil menghasilkan error yang disebut error pembulatan.

Sebagai contoh 1/6 = 0.166666666… tidak dapat dinyatakan secara tepat oleh komputer karena digit 6 panjangnya tidak terbatas. Komputer hanya mampu merepresentasikan sejumlah digit (atau bit dalam sistem biner) saja. Bilangan riil yang panjangnya melebihi jumlah digit (bit) yang dapat direpresentasikan oleh komputer dibulatkan ke bilangan terdekat. Misalnya sebuah komputer hanya dapat merepresentasikan bilangan riil dalam 6 digit angka berarti, maka representasi bilangan 1/6 = 0.1666666666… di dalam komputer 6-digit tersebut adalah 0.166667. Galat pembulatannya adalah 1/6 – 0.166667 = -0.000000333.

Contoh dalam sistem biner misalnya 1/10 = 0.000110011001100110011 00110011…2 direpresentasikan di dalam komputer dalam jumlah bit yang terbatas. Kebanyakan komputer digital mempunyai dua buah cara penyajian bilangan riil, yaitu bilangan titik-tetap (fixed point) dan bilangan titik-kambang (floatingpoint). Dalam format bilangan titik -tetap setiap bilangan disajikan dengan jumlah tempat desimal yang tetap, misalnya 62.358, 0.013, 1.000. Sedangkan dalam format bilangan titik-kambang setiap bilangan disajikan dengan jumlah digit berarti yang sudah tetap, misalnya 0.6238 103 0.1714 ^10-13 atau ditulis juga 0.6238E+03 0.1714E-13. Digit-digit berarti di dalam format bilangan titik-kambang disebut juga angka bena (significant figure).

**Galat Pemotongan ( truncation error )**

Galat pemotongan adalah galat yang ditimbulkan oleh pembatasan jumlah komputasi yang digunakan pada proses metode numerik. Banyak metode dalam metode numerik yang penurunan rumusnya menggunakan proses iterasi yang jumlahnya tak terhingga, sehingga untuk membatasi proses penghitungan, jumlah iterasi dibatasi sampai langkah ke n. Hasil penghitungan sampai langkah ke n akan menjadi hasil hampiran dan nilai penghitungan langkah n keatas akan menjadi galat pemotongan. dalam hal ini galat pemotongan kan menjadi sangat kecil sekali jika nilai n di perbesar. Konsekuensinya tentu saja jumlah proses penghitungannya akan semakin banyak.

**3. Cara atau penyebab lainnya antara lain :**
**1. Round off errors**
  Kesalahan ini biasanya akibat proses pembulatan dalam perhitungan. Secara umum, proses pembulatan ada 2 aturan yaitu :
-          Jika digit yang dibulatkan kurang dari 5, maka tidak terjadi pembulatan.
-          Sebaliknya, jika lebih dari 5, maka terjadi pembulatan yaitu dengan menambah satu.

**2. Truncation errors**
  Kesalahan pemotongan biasanya terjadi karena pembuangan suku yang berderajat tinggi. Sebagai contoh untuk menghitung nilai cosinus dapat menggunakan deret Taylor
  
**3.Range errors**
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
                
**Referensi**

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

**METODE BISEKSI (BISECTION METHOD)**
 
Dinamakan metode biseksi (Bi Section) didasarkan atas teknis metode ini adalah “belah dua”. Metode Biseksi diformulasikan berdasarkan Teorema konsep dasar kalkulus : nilai antara dan deret Taylor Teorema 1.1 yang menyatakan bahwa : “Bila fungsi f(x) kontinu dalam selang/interval (a,b) dengan f(a) dan f(b) berlawanan tanda, maka f(α)=0 untuk suatu bilangan α sedemikian hingga a<α<b.

Dinamakan metode biseksi (Bi Section) didasarkan atas teknis metode ini adalah “belah dua”. Metode Biseksi diformulasikan berdasarkan Teorema konsep dasar kalkulus : nilai antara dan deret Taylor Teorema 1.1 yang menyatakan bahwa : “Bila fungsi f(x) kontinu dalam selang/interval (a,b) dengan f(a) dan f(b) berlawanan tanda, maka f(α)=0 untuk suatu bilangan α sedemikian hingga a<α<b.

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/3.png" width="69" height="42">

Bila <img class="alignnone size-full wp-image-77" alt="rumus 4b" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b1.png" width="53" height="18"> atau <img class="alignnone size-full wp-image-78" alt="rumus 4b2" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b2.png" width="32" height="17"> mendekati nilai 0 untuk suatu nilai toleransi yang diberikan maka x0 adalah nilai akar dari <img class="alignnone size-full wp-image-77" alt="rumus 4b" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b1.png" width="53" height="18">

Sebaliknya bila <img class="alignnone size-full wp-image-77" alt="rumus 4b" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b1.png" width="53" height="18"> tidak sama dengan 0 atau <img class="alignnone size-full wp-image-78" alt="rumus 4b2" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-4b2.png" width="32" height="17"> mendekati nilai 0 tetapi **tidak memenuhi suatu nilai toleransi yang diberikan**, maka berdasarkan Teorema 1.1 ada dua kemungkinan yakni nilai akar berada di antara a dan x0 atau nilai akar berada di antara x0 dan b. Dari salah satu kemungkinan ini, metode Biseksi kembali akan digunakan. Secara geometris, metode Biseksi yang dikemukan di atas diilustrasikan melalui gambar grafik berikut ini.

<img class="alignnone size-full wp-image-79" alt="rumus 5" src="http://blog.ub.ac.id/aldipradana/files/2013/09/rumus-5.png" width="396" height="174">

**METODE NEWTON**

Metode NEWTON didasarkan pada aproksimasi linear fungsi dan menggunakan prinsip kemiringan (Tangen) kurvanya.

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/12.png" width="497" height="241">

Kalkulasi dengan metode Newton diawali dengan yang tidak terlalu jauh dari sebuah akar, bergerak sepanjang garis linear (kemiringan atau tangen garis) ke perpotongannya di sumbu-x, dan mengambilnya sebagai titik aproksimasi untuk yang berikutnya. Perlakuan ini diteruskan hingga nilai-nilai x dirasakan sukses cukup dekat ke fungsi bernilai nol. Skema kalkulasinya mengikuti segitiga yang dibangun dengan sudut inklinasi dari kemiringan garis pada kurva di yaitu

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/13.png" width="393" height="57">

Aproksimasi berikutnya diteruskan dengan menghitung  x2 dengan skema yang sama dimana nilai  x3 digantikan  oleh x1. Secara umum metode Newton dirumuskan oleh skema berikut ini:

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/14.png" width="235" height="54">

**METODE POSISI SALAH (REGULA FALSI)**

Untuk menghitung nilai akar dari f(x)=0 dapat digunakan metode Posisi salah/Regulasi Falsi. Aturan dari metode ini diterangkan secara geometri dalam Gambar. Dalam gambar yang dimaksud, sketsa grafik dari kurva dinyatakan oleh persamaan y=f(x). Akar dari f(x)=0 yang dicari dinyatakan oleh koordinat x dari titik P yang merupakan perpotongan dari kurva y=f(x) dengan sumbu x. Untuk menggunakan aturan RF, diperlukan dua titik A(α,f (α)) dan B0 (x0,f(x0)).

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/211-300x148.png" width="467" height="174">

Proses selanjutnya adalah menghitung nilai x melalui persamaan garis AB0 yang memotong sumbu x di titik P1. Setelah itu dengan menggunakan koordinat titik P1 yakni (x1, 0) dapat ditentukan titik B1 dengan koordinat (x1, f(x1)) . Dengan demikian garis AB1 akan memotong sumbu x di titik P2 dengan koordinat (x2, 0). Demikian proses ini terus dilakukan hingga diperoleh kondisi Pn sangat dekat ke P yakni |Pn-P|<toleransi. Dari proses pencapaian nilai akar di titik P, dihasilkan barisan nilai-nilai x0,x1,x2,……,xn yang diharapkan akan konvergen ke absis x pada titik P, yaitu akar yang dicari dari persamaan yang diberikan.

Persamaan garis AB0 adalah :

<img class="alignnone" alt="" src="http://blog.ub.ac.id/musthafaendybasranto/files/2013/09/22.png" width="244" height="61">

Karena persamaan garis AB0 melalui titik P1(x1,0), maka diperoleh :

<img class="alignnone size-medium wp-image-82" alt="a" src="http://blog.ub.ac.id/aldipradana/files/2013/09/a-300x145.jpg" width="243" height="116">

Sumber :

[http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&cad=rja&ved=0CEUQFjAD&url=http%3A%2F%2F
matematikaindo.
files.wordpress.com%2F2010%2F04%2Fmetode-
numerik-buku-ajar-unila.pdf&ei=cIA7UrhFi9KtB4-CgcgK&usg=AFQjCNGex_IcMyn18Q6O-ZmCHZ_E3638CQ&sig2=1C2ZFEGmJtwWu
zqhyEtFnQ](http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&cad=rja&ved=0CEUQFjAD&url=http%3A%2F%2Fmatematikaindo.files.wordpress.com%2F2010%2F04%2Fmetode-numerik-buku-ajar-unila.pdf&ei=cIA7UrhFi9KtB4-CgcgK&usg=AFQjCNGex_IcMyn18Q6O-ZmCHZ_E3638CQ&sig2=1C2ZFEGmJtwWuzqhyEtFnQ)

# ~ Numerical Differentiation

 Persamaan diferensial merupakan persoalan matematis yang sering dijumpai dalam bidang teknik lingkungan. Sering kali suatu persamaan diferensial tidak dapat diselesaikan secara analitik sehingga diperlukan metode numerik untuk menyelesaikannya.
 
 ~ Jenis-jenis Metode Diferensiasi : 
 Metode Euler
 Metode Runge-Kutta
 Metode Taylor
 Metode Adam
 Metode Milne
 Metode Adam-Moulton
 Metode modifikasi Euler

Differensial banyak digunakan dalam perhitungan kalkulus untuk keperluan perhitungan geometrik dan perubahan – perubahan nilai persatuan waktu atau jarak.

Differensial  merupakan perbandingan perubahan tinggi dan perubahan jarak yang secara kalkulus dapat didefinisikan sebagai berikut :

<img data-attachment-id="263" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/differential/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/differential.jpg" data-orig-size="134,68" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="differential" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/differential.jpg?w=134" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/differential.jpg?w=134" class=" size-full wp-image-263 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/differential.jpg?w=748" alt="differential">

Diferensiasi numerik digunakan dalam penentuan nilai aproksimasi dari turunan fungsi  f pada titik tertentu. Hubungan antara nilai fungsi dan perubahan fungsi untuk setiap titiknya didefinisikan sebagai berikut :

<img data-attachment-id="264" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/fungsi/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/fungsi.jpg" data-orig-size="241,64" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="fungsi" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/fungsi.jpg?w=241" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/fungsi.jpg?w=241" class=" size-full wp-image-264 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/fungsi.jpg?w=748" alt="fungsi.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/fungsi.jpg 241w, https://d4anm2017a.files.wordpress.com/2017/10/fungsi.jpg?w=150 150w" sizes="(max-width: 241px) 100vw, 241px">

dan f'(x) didefinisikan dengan :

<img data-attachment-id="265" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/dif/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/dif.jpg" data-orig-size="310,127" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="dif" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/dif.jpg?w=300" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/dif.jpg?w=310" class=" size-full wp-image-265 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/dif.jpg?w=748" alt="dif.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/dif.jpg 310w, https://d4anm2017a.files.wordpress.com/2017/10/dif.jpg?w=150 150w, https://d4anm2017a.files.wordpress.com/2017/10/dif.jpg?w=300 300w" sizes="(max-width: 310px) 100vw, 310px">

Terdapat 3 jenis diferensiasi dalam metode numerik yaitu :

**Metode Selisih Maju** 

Metode selisih maju merupakan metode yang mengadopsi secara langsung definisi differensial

<img data-attachment-id="266" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/selisih_maju/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju.jpg" data-orig-size="248,76" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="selisih_maju" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju.jpg?w=248" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju.jpg?w=248" class=" size-full wp-image-266 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju.jpg?w=748" alt="selisih_maju.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju.jpg 248w, https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju.jpg?w=150 150w" sizes="(max-width: 248px) 100vw, 248px">

Sehingga error yang dihasilkan

<img data-attachment-id="267" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/error_sel_maju/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_maju.jpg" data-orig-size="200,74" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="error_sel_maju" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_maju.jpg?w=200" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_maju.jpg?w=200" class=" size-full wp-image-267 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_maju.jpg?w=748" alt="error_sel_maju.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_maju.jpg 200w, https://d4anm2017a.files.wordpress.com/2017/10/error_sel_maju.jpg?w=150 150w" sizes="(max-width: 200px) 100vw, 200px">

Contoh :

Hitung differensial dari fungsi berikut

<img data-attachment-id="270" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/soal_sel_maju/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_maju-e1508385609994.jpg" data-orig-size="202,50" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="soal_sel_maju" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_maju-e1508385609994.jpg?w=202" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_maju-e1508385609994.jpg?w=202" class="alignnone size-full wp-image-270" src="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_maju-e1508385609994.jpg?w=748" alt="soal_sel_maju" srcset="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_maju-e1508385609994.jpg 202w, https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_maju-e1508385609994.jpg?w=150 150w" sizes="(max-width: 202px) 100vw, 202px">

dari range x = [0,1] dengan h = 0,05

<img data-attachment-id="271" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/ans_sel_maju/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/ans_sel_maju-e1508385721608.jpg" data-orig-size="483,320" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="ans_sel_maju" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/ans_sel_maju-e1508385721608.jpg?w=300" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/ans_sel_maju-e1508385721608.jpg?w=483" class=" size-full wp-image-271 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/ans_sel_maju-e1508385721608.jpg?w=748" alt="ans_sel_maju.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/ans_sel_maju-e1508385721608.jpg 483w, https://d4anm2017a.files.wordpress.com/2017/10/ans_sel_maju-e1508385721608.jpg?w=150 150w, https://d4anm2017a.files.wordpress.com/2017/10/ans_sel_maju-e1508385721608.jpg?w=300 300w" sizes="(max-width: 483px) 100vw, 483px">

**Metode Selisih Mundur**

Metode selisih mundur merupakan kebalikan dari metode selisih maju

Sehingga dapat didefinisikan dalam rumus berikut :

<img data-attachment-id="268" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/selisih_mundur/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_mundur.jpg" data-orig-size="303,193" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="selisih_mundur" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_mundur.jpg?w=300" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_mundur.jpg?w=303" class=" size-full wp-image-268 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/selisih_mundur.jpg?w=748" alt="selisih_mundur.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/selisih_mundur.jpg 303w, https://d4anm2017a.files.wordpress.com/2017/10/selisih_mundur.jpg?w=150 150w, https://d4anm2017a.files.wordpress.com/2017/10/selisih_mundur.jpg?w=300 300w" sizes="(max-width: 303px) 100vw, 303px">

**Metode Selisih tengah**

Metode selisih tengah merupakan metode pengambilan perubahan dari dua titik sekitar dari titik yang diukur.Perhatikan selisih maju pada titik x-h adalah :

<img data-attachment-id="272" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/selisih_maju_xh/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_xh.jpg" data-orig-size="218,80" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="selisih_maju_xh" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_xh.jpg?w=218" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_xh.jpg?w=218" class=" size-full wp-image-272 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_xh.jpg?w=748" alt="selisih_maju_xh.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_xh.jpg 218w, https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_xh.jpg?w=150 150w" sizes="(max-width: 218px) 100vw, 218px">

Dan selisih maju pada titik x adalah :

<img data-attachment-id="273" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/selisih_maju_x/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_x.jpg" data-orig-size="209,79" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="selisih_maju_x" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_x.jpg?w=209" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_x.jpg?w=209" class=" size-full wp-image-273 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_x.jpg?w=748" alt="selisih_maju_x.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_x.jpg 209w, https://d4anm2017a.files.wordpress.com/2017/10/selisih_maju_x.jpg?w=150 150w" sizes="(max-width: 209px) 100vw, 209px">

Metode selisih tengahan merupakan rata-rata dari dua selisih maju pada titik x-h dan titik x :

<img data-attachment-id="274" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/selisih_tengah/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_tengah.jpg" data-orig-size="275,139" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="selisih_tengah" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_tengah.jpg?w=275" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/selisih_tengah.jpg?w=275" class=" size-full wp-image-274 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/selisih_tengah.jpg?w=748" alt="selisih_tengah" srcset="https://d4anm2017a.files.wordpress.com/2017/10/selisih_tengah.jpg 275w, https://d4anm2017a.files.wordpress.com/2017/10/selisih_tengah.jpg?w=150 150w" sizes="(max-width: 275px) 100vw, 275px">

Sehingga error yang dihasilkan

<img data-attachment-id="277" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/error_sel_tengah/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_tengah.jpg" data-orig-size="195,82" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="error_sel_tengah" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_tengah.jpg?w=195" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_tengah.jpg?w=195" class="aligncenter size-full wp-image-277" src="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_tengah.jpg?w=748" alt="error_sel_tengah.jpg" srcset="https://d4anm2017a.files.wordpress.com/2017/10/error_sel_tengah.jpg 195w, https://d4anm2017a.files.wordpress.com/2017/10/error_sel_tengah.jpg?w=150 150w" sizes="(max-width: 195px) 100vw, 195px">

<img loading="lazy" data-attachment-id="278" data-permalink="https://d4anm2017a.wordpress.com/2017/10/19/diferensiasi-numerik/soal_sel_tengah/" data-orig-file="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg" data-orig-size="847,580" data-comments-opened="1" data-image-meta="{&quot;aperture&quot;:&quot;0&quot;,&quot;credit&quot;:&quot;&quot;,&quot;camera&quot;:&quot;&quot;,&quot;caption&quot;:&quot;&quot;,&quot;created_timestamp&quot;:&quot;0&quot;,&quot;copyright&quot;:&quot;&quot;,&quot;focal_length&quot;:&quot;0&quot;,&quot;iso&quot;:&quot;0&quot;,&quot;shutter_speed&quot;:&quot;0&quot;,&quot;title&quot;:&quot;&quot;,&quot;orientation&quot;:&quot;0&quot;}" data-image-title="soal_sel_tengah" data-image-description="" data-medium-file="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg?w=300" data-large-file="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg?w=748" class="  wp-image-278 aligncenter" src="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg?w=529&amp;h=363" alt="soal_sel_tengah.jpg" width="529" height="363" srcset="https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg?w=529&amp;h=363 529w, https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg?w=150&amp;h=103 150w, https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg?w=300&amp;h=205 300w, https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg?w=768&amp;h=526 768w, https://d4anm2017a.files.wordpress.com/2017/10/soal_sel_tengah.jpg 847w" sizes="(max-width: 529px) 100vw, 529px">

# ~ Numerical Integration

Metode integrasi numerik adalah suatu cara untuk menghitung aproksimasi luas daerah di bawah fungsi yang dimaksud pada selang yang diberikan. Berikut ini adalah beberapa metode integrasi numerik yang lazim digunakan :

**Metode Euler Eksplisit**
merupakan metode integrasi yang paling mudah 
<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/ce179f1ee061e052f254a491a9dbe7512e364d41" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:39.261ex; height:2.843ex;" alt="{\displaystyle {\dot {x}}_{k-1}=Ax_{k-1}+Bu_{k-1}=f(x_{k-1},u_{k-1})}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/17f33d249eefbffdc3d2919e7821782622d925cb" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.671ex; width:18.734ex; height:2.509ex;" alt="{\displaystyle x_{k}=x_{k-1}+h{\dot {x}}_{k-1}}">

**Metode Euler Implisit**
<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/5dc9446a2dc4b9560f842c945741e093ced41643" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:30.859ex; height:2.843ex;" alt="{\displaystyle {\dot {x}}_{k-1}=Ax_{k}+Bu_{k}=f(x_{k},u_{k})}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/f71a09b7421f85ef000b88a1afaff9172b9dd574" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.671ex; width:16.634ex; height:2.509ex;" alt="{\displaystyle x_{k}=x_{k-1}+h{\dot {x}}_{k}}">

Pada metode integrasi implisit nilai aktual <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/6d2b88c64c76a03611549fb9b4cf4ed060b56002" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.671ex; width:2.418ex; height:2.009ex;" alt="{\displaystyle x_{k}}"> juga digunakan sebagai umpan balik. Umpan balik ini dapat menyebabkan terjadinya lingkaran aljabar. Untuk menghindarinya maka bentuk persamaan diubah menjadi seperti ini

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/5e3797b09546e064dff2182315ffc4915e583573" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:32.959ex; height:2.843ex;" alt="{\displaystyle {\dot {x}}_{k}=Ax_{k-1}+Bu_{k}=f(x_{k-1},u_{k})}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/5c8d44a52caca3e6411978102f76493a5663e077" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:27.083ex; height:3.176ex;" alt="{\displaystyle x_{k}=x_{k-1}+h[I-hJ]^{-1}{\dot {x}}_{k}}"> 

J adalah matriks Jacobi. Pada sistem linear dan invarian terhadap waktu, maka matriks J = A

contoh program lain di dalam metode euler 
```markdown
euler <- function(f, x0, y0, h, n){
  x <- x0
  y <- y0
  
  for(i in 1:n){
    y0 <- y0 + h*f(x0, y0)
    x0 <- x0 + h
    x <- c(x,x0)
    y <- c(y, y0)
  }
  
  return(data.frame(x=x, y=y))
}
```

```markdown
# metode numerik
f1 <- function(x,y){y/(2*x+1)}
num <- euler(f1, x0=0, y0=1, h=0.05, n=100)

# metode analitik
f2 <- function(x){sqrt(2*x+1)}
x0 <- 0
y0 <- 1
x <- x0
y <- y0

for(i in 1:100){
  y0 <- f2(x0+0.05)
  x0 <- x0+0.05
  x <- c(x, x0)
  y <- c(y, y0)
}
true <- data.frame(x=x, y=y)
```

**Metode Heun**
Algoritma integrasi Heun memerlukan dua masukan yaitu <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/9837644700489d04d977da272524cd5fda36f3d7" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.671ex; width:2.418ex; height:2.009ex;" alt="{\displaystyle u_{k}}"> dan <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/3ab4fbdd181a4697fd87bb4b73926d99b1bc8d59" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.671ex; width:4.519ex; height:2.009ex;" alt="{\displaystyle u_{k-1}}"> 

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/ce179f1ee061e052f254a491a9dbe7512e364d41" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:39.261ex; height:2.843ex;" alt="{\displaystyle {\dot {x}}_{k-1}=Ax_{k-1}+Bu_{k-1}=f(x_{k-1},u_{k-1})}">

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/2daa0dd975f2173e23a3d1d92e890e6dd0b7f047" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.005ex; width:18.734ex; height:3.176ex;" alt="{\displaystyle x_{k}^{p}=x_{k-1}+h{\dot {x}}_{k-1}}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/a820e9bf3ffa92a3753fafa5880a3160335cfa27" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.005ex; width:14.476ex; height:3.176ex;" alt="{\displaystyle {\dot {x}}_{k}^{p}=f(x_{k}^{p},u_{k})}">

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/7fa4aab434af5944e0dd7aefe53bc6f90d80517e" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:26.638ex; height:5.343ex;" alt="{\displaystyle x_{k}=x_{k-1}+{h \over 2}({\dot {x}}_{k-1}+{\dot {x}}_{k}^{p})}">

contoh program metode ini

```markdown
heun <- function(f, x0, y0, h, n, iter=1){
  x <- x0
  y <- y0
  
  for(i in 1:n){
    ypred0 <- f(x0,y0)
    ypred1 <- y0 + h*ypred0
    ypred2 <- f(x0+h,ypred1)
    ykor <- y0 + h*(ypred0+ypred2)/2
    if(iter!=1){
      for(i in 1:iter){
        ykor <- y0 + h*(ypred0+f(x0+h,ykor))/2
      }
    }
    y0 <- ykor
    x0 <- x0 + h
    x <- c(x, x0)
    y <- c(y, y0)
  }
  
  return(data.frame(x=x,y=y))
}
```

**Metode Runge-Kutta**

Merupakan integrator dengan empat masukan.

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/ce179f1ee061e052f254a491a9dbe7512e364d41" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:39.261ex; height:2.843ex;" alt="{\displaystyle {\dot {x}}_{k-1}=Ax_{k-1}+Bu_{k-1}=f(x_{k-1},u_{k-1})}">

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/d92d90efaa36b7b54177554d63cc1f581996b7e9" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:22.95ex; height:5.343ex;" alt="{\displaystyle x_{k-0.5}^{p1}=x_{k-1}+{h \over 2}{\dot {x}}_{k-1}}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/280f617f3e84dd807d202a2226a4daf2c5018e22" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.171ex; width:24.615ex; height:3.676ex;" alt="{\displaystyle {\dot {x}}_{k-0.5}^{p1}=f(x_{k-0.5}^{p1},u_{k-0.5})}">

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/c380357ae13093517fd6f97d652ea5458224e044" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:24.23ex; height:5.343ex;" alt="{\displaystyle x_{k-0.5}^{p2}=x_{k-1}+{h \over 2}{\dot {x}}_{k-0.5}^{p1}}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/382b9580ac0ed989755a6bf13a98f6173f47e8ed" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.171ex; width:24.615ex; height:3.676ex;" alt="{\displaystyle {\dot {x}}_{k-0.5}^{p2}=f(x_{k-0.5}^{p2},u_{k-0.5})}">

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/d97472d30e042c7eddbc6047301a226771aa8847" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.171ex; width:20.806ex; height:3.676ex;" alt="{\displaystyle x_{k}^{p3}=x_{k-1}+h{\dot {x}}_{k-0.5}^{p2}}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/e6fdc1ddc1e929c5e1ae7a70f2f49cdfe9c69379" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.005ex; width:16.06ex; height:3.509ex;" alt="{\displaystyle {\dot {x}}_{k}^{p3}=f(x_{k}^{p3},u_{k})}">

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/0389ac3fd60b7f9113005ab609caa69cf1bad324" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:47.033ex; height:5.343ex;" alt="{\displaystyle x_{k}=x_{k-1}+{h \over 6}({\dot {x}}_{k-1}+2{\dot {x}}_{k-0.5}^{p1}+2{\dot {x}}_{k-0.5}^{p2}+{\dot {x}}_{k}^{p3})}"> 

**Metode Trapesium (Trapez)** 

merupakan nilai tengah dari metode Euler eksplisit dan metode Euler implisit.

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/ce179f1ee061e052f254a491a9dbe7512e364d41" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:39.261ex; height:2.843ex;" alt="{\displaystyle {\dot {x}}_{k-1}=Ax_{k-1}+Bu_{k-1}=f(x_{k-1},u_{k-1})}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/3ece1f24cd7f13bf2488ec2e5cf1fdefc0f772a3" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:28.758ex; height:2.843ex;" alt="{\displaystyle {\dot {x}}_{k}=Ax_{k}+Bu_{k}=f(x_{k},u_{k})}"> 

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/c64069d0a7d341e11390cac1178b5856af6d7a09" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:26.638ex; height:5.343ex;" alt="{\displaystyle x_{k}=x_{k-1}+{h \over 2}({\dot {x}}_{k}+{\dot {x}}_{k+1})}">

Sama halnya dengan metode Euler implisit, metode ini dapat menyebabkan lingkaran aljabar. Oleh karena itu, bentuk persamaan ini diubah menjadi seperti ini

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/fea32d9c7dea988bb29126e2a375e17417275fe5" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:50.618ex; height:5.176ex;" alt="{\displaystyle {\dot {x}}_{k-1}=Ax_{k-1}+{B \over 2}(u_{k-1}+u_{k})=f(x_{k-1},u_{k-1},u_{k})}"> <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/d8ecf44f47bdee29bc63d9a1b9da8592a82f690e" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:27.919ex; height:5.343ex;" alt="{\displaystyle x_{k}=x_{k-1}+h[I-{h \over 2}J]^{-1}{\dot {x}}_{k}}"> 

**Metode Newton–Cotes**

<table class="wikitable" style="margin:1em auto 1em auto; background:white">

<tbody><tr>
<th>No.</th>
<th>Nama Aturan</th>
<th>Rumus</th>
<th>Estimasi Kesalahan
</th></tr>
<tr align="center">
<td>1</td>
<td><a href="/wiki/Trapesium_(geometri)" title="Trapesium (geometri)">Trapesium</a></td>
<td><span class="mwe-math-element"><span class="mwe-math-mathml-inline mwe-math-mathml-a11y" style="display: none;"><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="{\displaystyle {\frac {b-a}{2}}(f_{0}+f_{1})}">
  <semantics>
    <mrow class="MJX-TeXAtom-ORD">
      <mstyle displaystyle="true" scriptlevel="0">
        <mrow class="MJX-TeXAtom-ORD">
          <mfrac>
            <mrow>
              <mi>b</mi>
              <mo>−<!-- − --></mo>
              <mi>a</mi>
            </mrow>
            <mn>2</mn>
          </mfrac>
        </mrow>
        <mo stretchy="false">(</mo>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>0</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>1</mn>
          </mrow>
        </msub>
        <mo stretchy="false">)</mo>
      </mstyle>
    </mrow>
    <annotation encoding="application/x-tex">{\displaystyle {\frac {b-a}{2}}(f_{0}+f_{1})}</annotation>
  </semantics>
</math></span><img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/898874be93dd6abb959fcbef070b3ddcb08c96f9" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:14.941ex; height:5.343ex;" alt="{\displaystyle {\frac {b-a}{2}}(f_{0}+f_{1})}"></span></td>
<td><span class="mwe-math-element"><span class="mwe-math-mathml-inline mwe-math-mathml-a11y" style="display: none;"><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="{\displaystyle -{\frac {(b-a)^{3}}{12}}\,f^{(2)}(\xi )}">
  <semantics>
    <mrow class="MJX-TeXAtom-ORD">
      <mstyle displaystyle="true" scriptlevel="0">
        <mo>−<!-- − --></mo>
        <mrow class="MJX-TeXAtom-ORD">
          <mfrac>
            <mrow>
              <mo stretchy="false">(</mo>
              <mi>b</mi>
              <mo>−<!-- − --></mo>
              <mi>a</mi>
              <msup>
                <mo stretchy="false">)</mo>
                <mrow class="MJX-TeXAtom-ORD">
                  <mn>3</mn>
                </mrow>
              </msup>
            </mrow>
            <mn>12</mn>
          </mfrac>
        </mrow>
        <mspace width="thinmathspace"></mspace>
        <msup>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mo stretchy="false">(</mo>
            <mn>2</mn>
            <mo stretchy="false">)</mo>
          </mrow>
        </msup>
        <mo stretchy="false">(</mo>
        <mi>ξ<!-- ξ --></mi>
        <mo stretchy="false">)</mo>
      </mstyle>
    </mrow>
    <annotation encoding="application/x-tex">{\displaystyle -{\frac {(b-a)^{3}}{12}}\,f^{(2)}(\xi )}</annotation>
  </semantics>
</math></span><img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/832a38281b02d44403f66937c749653209eefe07" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:17.456ex; height:5.843ex;" alt="{\displaystyle -{\frac {(b-a)^{3}}{12}}\,f^{(2)}(\xi )}"></span>
</td></tr>
<tr align="center">
<td>2</td>
<td><a href="/wiki/Kaidah_Simpson" title="Kaidah Simpson">Simpson 1/3</a></td>
<td><span class="mwe-math-element"><span class="mwe-math-mathml-inline mwe-math-mathml-a11y" style="display: none;"><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="{\displaystyle {\frac {b-a}{3}}(f_{0}+4f_{1}+f_{2})}">
  <semantics>
    <mrow class="MJX-TeXAtom-ORD">
      <mstyle displaystyle="true" scriptlevel="0">
        <mrow class="MJX-TeXAtom-ORD">
          <mfrac>
            <mrow>
              <mi>b</mi>
              <mo>−<!-- − --></mo>
              <mi>a</mi>
            </mrow>
            <mn>3</mn>
          </mfrac>
        </mrow>
        <mo stretchy="false">(</mo>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>0</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <mn>4</mn>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>1</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>2</mn>
          </mrow>
        </msub>
        <mo stretchy="false">)</mo>
      </mstyle>
    </mrow>
    <annotation encoding="application/x-tex">{\displaystyle {\frac {b-a}{3}}(f_{0}+4f_{1}+f_{2})}</annotation>
  </semantics>
</math></span><img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/f0e59bd2336d992fca79d6bc72be3cc75f59a350" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:21.137ex; height:5.343ex;" alt="{\displaystyle {\frac {b-a}{3}}(f_{0}+4f_{1}+f_{2})}"></span></td>
<td><span class="mwe-math-element"><span class="mwe-math-mathml-inline mwe-math-mathml-a11y" style="display: none;"><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="{\displaystyle -{\frac {(b-a)^{5}}{90}}\,f^{(4)}(\xi )}">
  <semantics>
    <mrow class="MJX-TeXAtom-ORD">
      <mstyle displaystyle="true" scriptlevel="0">
        <mo>−<!-- − --></mo>
        <mrow class="MJX-TeXAtom-ORD">
          <mfrac>
            <mrow>
              <mo stretchy="false">(</mo>
              <mi>b</mi>
              <mo>−<!-- − --></mo>
              <mi>a</mi>
              <msup>
                <mo stretchy="false">)</mo>
                <mrow class="MJX-TeXAtom-ORD">
                  <mn>5</mn>
                </mrow>
              </msup>
            </mrow>
            <mn>90</mn>
          </mfrac>
        </mrow>
        <mspace width="thinmathspace"></mspace>
        <msup>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mo stretchy="false">(</mo>
            <mn>4</mn>
            <mo stretchy="false">)</mo>
          </mrow>
        </msup>
        <mo stretchy="false">(</mo>
        <mi>ξ<!-- ξ --></mi>
        <mo stretchy="false">)</mo>
      </mstyle>
    </mrow>
    <annotation encoding="application/x-tex">{\displaystyle -{\frac {(b-a)^{5}}{90}}\,f^{(4)}(\xi )}</annotation>
  </semantics>
</math></span><img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/e6c3700425a20e7b2f6975207f689fd8cc139d09" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:17.456ex; height:5.843ex;" alt="{\displaystyle -{\frac {(b-a)^{5}}{90}}\,f^{(4)}(\xi )}"></span>
</td></tr>
<tr align="center">
<td>3</td>
<td><a href="/wiki/Kaidah_Simpson" title="Kaidah Simpson">Simpson 3/8</a></td>
<td><span class="mwe-math-element"><span class="mwe-math-mathml-inline mwe-math-mathml-a11y" style="display: none;"><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="{\displaystyle {\frac {3(b-a)}{8}}(f_{0}+3f_{1}+3f_{2}+f_{3})}">
  <semantics>
    <mrow class="MJX-TeXAtom-ORD">
      <mstyle displaystyle="true" scriptlevel="0">
        <mrow class="MJX-TeXAtom-ORD">
          <mfrac>
            <mrow>
              <mn>3</mn>
              <mo stretchy="false">(</mo>
              <mi>b</mi>
              <mo>−<!-- − --></mo>
              <mi>a</mi>
              <mo stretchy="false">)</mo>
            </mrow>
            <mn>8</mn>
          </mfrac>
        </mrow>
        <mo stretchy="false">(</mo>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>0</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <mn>3</mn>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>1</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <mn>3</mn>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>2</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>3</mn>
          </mrow>
        </msub>
        <mo stretchy="false">)</mo>
      </mstyle>
    </mrow>
    <annotation encoding="application/x-tex">{\displaystyle {\frac {3(b-a)}{8}}(f_{0}+3f_{1}+3f_{2}+f_{3})}</annotation>
  </semantics>
</math></span><img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/7142aae98d40aa460f706fb825531ff99883d5c1" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:30.305ex; height:5.676ex;" alt="{\displaystyle {\frac {3(b-a)}{8}}(f_{0}+3f_{1}+3f_{2}+f_{3})}"></span></td>
<td><span class="mwe-math-element"><span class="mwe-math-mathml-inline mwe-math-mathml-a11y" style="display: none;"><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="{\displaystyle -{\frac {3(b-a)^{5}}{80}}\,f^{(4)}(\xi )}">
  <semantics>
    <mrow class="MJX-TeXAtom-ORD">
      <mstyle displaystyle="true" scriptlevel="0">
        <mo>−<!-- − --></mo>
        <mrow class="MJX-TeXAtom-ORD">
          <mfrac>
            <mrow>
              <mn>3</mn>
              <mo stretchy="false">(</mo>
              <mi>b</mi>
              <mo>−<!-- − --></mo>
              <mi>a</mi>
              <msup>
                <mo stretchy="false">)</mo>
                <mrow class="MJX-TeXAtom-ORD">
                  <mn>5</mn>
                </mrow>
              </msup>
            </mrow>
            <mn>80</mn>
          </mfrac>
        </mrow>
        <mspace width="thinmathspace"></mspace>
        <msup>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mo stretchy="false">(</mo>
            <mn>4</mn>
            <mo stretchy="false">)</mo>
          </mrow>
        </msup>
        <mo stretchy="false">(</mo>
        <mi>ξ<!-- ξ --></mi>
        <mo stretchy="false">)</mo>
      </mstyle>
    </mrow>
    <annotation encoding="application/x-tex">{\displaystyle -{\frac {3(b-a)^{5}}{80}}\,f^{(4)}(\xi )}</annotation>
  </semantics>
</math></span><img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/ed57b1cd5a2f3daf1500910afccd87ae7669db95" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.838ex; width:18.619ex; height:5.843ex;" alt="{\displaystyle -{\frac {3(b-a)^{5}}{80}}\,f^{(4)}(\xi )}"></span>
</td></tr>
<tr align="center">
<td>4</td>
<td><a href="/wiki/Boole" class="mw-redirect" title="Boole">Boole</a> atau <a href="/wiki/Bode" class="mw-redirect" title="Bode">Bode</a></td>
<td><span class="mwe-math-element"><span class="mwe-math-mathml-inline mwe-math-mathml-a11y" style="display: none;"><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="{\displaystyle {\frac {2(b-a)}{45}}(7f_{0}+32f_{1}+12f_{2}+32f_{3}+7f_{4})}">
  <semantics>
    <mrow class="MJX-TeXAtom-ORD">
      <mstyle displaystyle="true" scriptlevel="0">
        <mrow class="MJX-TeXAtom-ORD">
          <mfrac>
            <mrow>
              <mn>2</mn>
              <mo stretchy="false">(</mo>
              <mi>b</mi>
              <mo>−<!-- − --></mo>
              <mi>a</mi>
              <mo stretchy="false">)</mo>
            </mrow>
            <mn>45</mn>
          </mfrac>
        </mrow>
        <mo stretchy="false">(</mo>
        <mn>7</mn>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>0</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <mn>32</mn>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>1</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <mn>12</mn>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>2</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <mn>32</mn>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>3</mn>
          </mrow>
        </msub>
        <mo>+</mo>
        <mn>7</mn>
        <msub>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mn>4</mn>
          </mrow>
        </msub>
        <mo stretchy="false">)</mo>
      </mstyle>
    </mrow>
    <annotation encoding="application/x-tex">{\displaystyle {\frac {2(b-a)}{45}}(7f_{0}+32f_{1}+12f_{2}+32f_{3}+7f_{4})}</annotation>
  </semantics>
</math></span><img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/0aac892bb48dad85a8538d2b1f4006bf0dff2cf9" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -2.005ex; width:42.314ex; height:5.843ex;" alt="{\displaystyle {\frac {2(b-a)}{45}}(7f_{0}+32f_{1}+12f_{2}+32f_{3}+7f_{4})}"></span></td>
<td><span class="mwe-math-element"><span class="mwe-math-mathml-inline mwe-math-mathml-a11y" style="display: none;"><math xmlns="http://www.w3.org/1998/Math/MathML" alttext="{\displaystyle -{\frac {8(b-a)^{7}}{945}}\,f^{(6)}(\xi )}">
  <semantics>
    <mrow class="MJX-TeXAtom-ORD">
      <mstyle displaystyle="true" scriptlevel="0">
        <mo>−<!-- − --></mo>
        <mrow class="MJX-TeXAtom-ORD">
          <mfrac>
            <mrow>
              <mn>8</mn>
              <mo stretchy="false">(</mo>
              <mi>b</mi>
              <mo>−<!-- − --></mo>
              <mi>a</mi>
              <msup>
                <mo stretchy="false">)</mo>
                <mrow class="MJX-TeXAtom-ORD">
                  <mn>7</mn>
                </mrow>
              </msup>
            </mrow>
            <mn>945</mn>
          </mfrac>
        </mrow>
        <mspace width="thinmathspace"></mspace>
        <msup>
          <mi>f</mi>
          <mrow class="MJX-TeXAtom-ORD">
            <mo stretchy="false">(</mo>
            <mn>6</mn>
            <mo stretchy="false">)</mo>
          </mrow>
        </msup>
        <mo stretchy="false">(</mo>
        <mi>ξ<!-- ξ --></mi>
        <mo stretchy="false">)</mo>
      </mstyle>
    </mrow>
    <annotation encoding="application/x-tex">{\displaystyle -{\frac {8(b-a)^{7}}{945}}\,f^{(6)}(\xi )}</annotation>
  </semantics>
</math></span><img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/0299d1ce872031ad9c7295e7a769c100d4c8d57e" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -2.005ex; width:18.619ex; height:6.009ex;" alt="{\displaystyle -{\frac {8(b-a)^{7}}{945}}\,f^{(6)}(\xi )}"></span>
</td></tr></tbody></table>





# ~ Numerical Solution of System of Linear Algebraic Ecuations

Sistem Persamaan Aljabar Linier (SPAL) atau dikenal juga sebagai ‘Persamaan   Aljabar   Linier   Serempak’   banyak   sekali   dijumpai dalam   perhitungan-perhitungan   teknik   kimia   yang   melibatkan solusi  numeris.  Beberapa  metode  solusi  yang  melibatkan  solusi SPAL,  di  antaranya  dalah:  solusi  Sisten  Persamaan  Aljabar  Non-Linier (SPANL), solusi Persamaan Diferensial Biasa (PDB), solusi persamaan  Diferensial  Parsial  (PDP),  Regresi  Linier  dan  Non-Linier, dll.

Macam - macam Metode Numeriknya yaitu :

**Eliminasi Gauss** 

- Eliminasi Gauss adalah suatu metode untuk mengoperasikan nilai-nilai di dalam matriks sehingga menjadi matriks yang lebih sederhana lagi. 
- Dengan melakukan operasi baris sehingga matriks tersebut menjadi matriks yang baris. Ini dapat digunakan sebagai salah satu metode penyelesaian persamaan linear dengan menggunakan matriks. 
- Caranya dengan mengubah persamaan linear tersebut ke dalam matriks teraugmentasi dan mengoperasikannya. Setelah menjadi matriks baris, lakukan substitusi balik untuk mendapatkan nilai dari variabel-variabel tersebut.

**ciri ciri Eliminasi Gauss**

- Jika suatu baris tidak semua nol, maka bilangan pertama yang tidak nol adalah 1 (1 utama)
- Baris nol terletak paling bawah
- 1 utama baris berikutnya berada dikanan 1 utama baris diatasnya
- Dibawah 1 utama harus nol.

**cara Pengerjaan Eliminasi Gauss**
- Metode Eliminasi Gauss merupakan metode yang dikembangkan dari metode eliminasi, yaitu menghilangkan atau mengurangi jumlah variable sehingga dapat diperoleh nilai dari suatu variable bebas
- matrik diubah menjadi augmented matrik :

<img border="0" data-original-height="253" data-original-width="353" height="229" src="https://1.bp.blogspot.com/-A5WDpPa0Trc/Xju-c7lXN2I/AAAAAAAAadI/ylmX6aa9A0gOe81G1_uv2PHmOCHeozECQCNcBGAsYHQ/s320/Slide5.jpg" width="320">

 ubah matrik menjadi matrik segitiga atas atau segitiga bawah dengan menggunakan OBE (Operasi Baris Elementer).
 
 <img border="0" data-original-height="260" data-original-width="918" height="179" src="https://1.bp.blogspot.com/-qQyodaiso7A/Xju-ci5UGXI/AAAAAAAAadE/Df8RikWTfX06l7GWHbS_iR4HLa2k5UUegCNcBGAsYHQ/s640/Slide6.jpg" width="640">
 
 **Operasi Baris Elementer**
 Metode dasar untuk menyelesaikan Sistem Persamaan Linier adalah mengganti sistem yang ada dengan sistem yang baru yang mempunyai himp solusi yang sama dan lebih mudah untuk diselesaikan

Sistem yang baru diperoleh dengan serangkaian step yang menerapkan 3 tipe operasi. Operasi ini disebut Operasi Baris Elementer

1. kalikan persamaan dengan nilai konstan yang bukan nol.
2. menukar dua persamaan.
3. menambahkan kelipatan dengan persamaan lain.

<img border="0" data-original-height="720" data-original-width="960" height="480" src="https://1.bp.blogspot.com/-qpqLXRzNKrs/Xju_q4SLxWI/AAAAAAAAadc/T4RGWmdCZdgQXm8WkTw1yWGxZ2w823tKACNcBGAsYHQ/s640/Slide9.JPG" width="640">

<img border="0" data-original-height="720" data-original-width="960" height="480" src="https://1.bp.blogspot.com/-NiHfRC-481k/Xju_qqSXRyI/AAAAAAAAadY/zm5UYrwQRUc0SOhTTVHM7-LJL1crp4eFQCNcBGAsYHQ/s640/Slide10.JPG" width="640">

<img border="0" data-original-height="720" data-original-width="960" height="480" src="https://1.bp.blogspot.com/-R9x2uA24PYk/Xju_qwTOhTI/AAAAAAAAadg/vr93e2sQ3UE2dkJsRsT0SVSna6Aeqvl2ACNcBGAsYHQ/s640/Slide11.JPG" width="640">

# ~ Numerical Solution of Ordinary Differential Equations

Salah satu dari beberapa jenis dalam sub bab ini yaitu dengan metode euler/metode heunn

Metode Euler digunakan sebagai dasar untuk metode Heun. Metode Euler menggunakan garis tangen ke fungsi di awal interval sebagai perkiraan kemiringan fungsi selama interval, dengan asumsi bahwa jika ukuran langkah kecil, kesalahan akan menjadi kecil. Namun, bahkan ketika ukuran langkah sangat kecil digunakan, lebih dari sejumlah besar langkah kesalahan mulai menumpuk dan estimasi menyimpang dari nilai fungsional aktual.

Di mana kurva solusi cekung ke atas, garis singgungnya akan meremehkan koordinat vertikal dari titik berikutnya dan sebaliknya untuk solusi cekung ke bawah. Garis prediksi ideal akan mengenai kurva pada titik prediksi berikutnya. Pada kenyataannya, tidak ada cara untuk mengetahui apakah solusinya cekung atau cekung, dan karenanya jika titik prediksi selanjutnya akan melebih-lebihkan atau meremehkan nilai vertikalnya. Konkavitas kurva tidak dapat dijamin untuk tetap konsisten dan prediksi dapat melebih-lebihkan dan meremehkan pada titik yang berbeda dalam domain solusi. Metode Heun mengatasi masalah ini dengan mempertimbangkan interval yang direntang oleh segmen garis singgung secara keseluruhan. Mengambil contoh cekung, garis prediksi tangen kiri meremehkan kemiringan kurva untuk seluruh lebar interval dari titik saat ini ke titik prediksi berikutnya. Jika garis singgung pada titik ujung kanan dianggap (yang dapat diperkirakan dengan menggunakan Metode Euler), ia memiliki masalah yang berlawanan[3] Titik-titik di sepanjang garis singgung dari titik ujung kiri memiliki koordinat vertikal yang semuanya meremehkan orang-orang yang bersandar pada kurva solusi, termasuk titik ujung kanan interval yang dipertimbangkan. Solusinya adalah membuat kemiringan lebih besar dengan jumlah tertentu. Metode Heun mempertimbangkan garis singgung ke kurva solusi di kedua ujung interval, yang terlalu tinggi , dan yang meremehkan koordinat vertikal yang ideal. Garis prediksi harus dibangun berdasarkan kemiringan tangen titik ujung kanan saja, diperkirakan menggunakan Metode Euler. Jika kemiringan ini melewati titik ujung kiri dari interval, hasilnya jelas terlalu curam untuk digunakan sebagai garis prediksi ideal dan melebih-lebihkan titik ideal. Oleh karena itu, titik ideal terletak kira-kira di tengah-tengah antara perkiraan yang salah dan terlalu rendah, rata-rata dari dua lereng.

Metode Euler digunakan untuk memperkirakan secara kasar koordinat titik berikutnya dalam solusi, dan dengan pengetahuan ini, perkiraan awal diprediksi ulang atau diperbaiki.
Dengan asumsi kuantitas <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/35b1a080f69eda39aa8be775b628dcb3a0c1919b" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:6.607ex; height:2.843ex;" alt="{\displaystyle \textstyle f(x,y)}"> di sisi kanan persamaan dapat dianggap sebagai kemiringan dari solusi yang dicari di titik mana pun <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/534d6f5286e49bcbb0120b3529b6043ba3b84d80" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:5.328ex; height:2.843ex;" alt="{\displaystyle \textstyle (x,y)}">  ini dapat dikombinasikan dengan estimasi Euler dari titik berikutnya untuk memberikan kemiringan garis singgung di titik akhir kanan. Selanjutnya rata-rata kedua lereng digunakan untuk menemukan koordinat yang diperbaiki dari interval ujung kanan. 

<img alt="Heun's Method." src="//upload.wikimedia.org/wikipedia/commons/thumb/5/52/Heun%27s_Method_Diagram.jpg/220px-Heun%27s_Method_Diagram.jpg" decoding="async" width="220" height="165" class="thumbimage" srcset="//upload.wikimedia.org/wikipedia/commons/thumb/5/52/Heun%27s_Method_Diagram.jpg/330px-Heun%27s_Method_Diagram.jpg 1.5x, //upload.wikimedia.org/wikipedia/commons/thumb/5/52/Heun%27s_Method_Diagram.jpg/440px-Heun%27s_Method_Diagram.jpg 2x" data-file-width="1024" data-file-height="768">

Prosedur untuk menghitung solusi numerik untuk masalah nilai awal: 

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/c23fb210b5b081d30b297cb7e7ec54e9522e7cd7" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:37.48ex; height:3.009ex;" alt="{\displaystyle y'(t)=f(t,y(t)),\qquad \qquad y(t_{0})=y_{0},}">

dengan cara metode Heun, adalah pertama menghitung nilai antara <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/c23fb210b5b081d30b297cb7e7ec54e9522e7cd7" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.838ex; width:37.48ex; height:3.009ex;" alt="{\displaystyle y'(t)=f(t,y(t)),\qquad \qquad y(t_{0})=y_{0},}"> dan kemudian perkiraan akhir <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/b92fefd5a55d5006605e793464e0fd56f6e13a3d" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.671ex; width:4.039ex; height:2.009ex;" alt="{\displaystyle y_{i+1}}"> pada titik integrasi berikutnya.

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/28dd4e6a8b0d676be1ad847893ee73608877e696" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -1.005ex; width:21.119ex; height:3.009ex;" alt="{\displaystyle {\tilde {y}}_{i+1}=y_{i}+hf(t_{i},y_{i})}">

dimana <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/b26be3e694314bc90c3215047e4a2010c6ee184a" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.338ex; width:1.339ex; height:2.176ex;" alt="{\displaystyle h}"> adalah ukuran langkah dan <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/8326d2c5936361393522689d49de51e49dfa7c61" class="mwe-math-fallback-image-inline" aria-hidden="true" style="vertical-align: -0.671ex; width:12.657ex; height:2.509ex;" alt="{\displaystyle t_{i+1}=t_{i}+h}"> 
















