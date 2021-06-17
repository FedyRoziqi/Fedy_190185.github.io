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


