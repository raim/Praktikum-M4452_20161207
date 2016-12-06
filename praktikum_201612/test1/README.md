# prepare layout file

1. export to .csv with tab as delimiter
2. replace underscore between I and aTc fields:
sed -i 's/I_/I;/g' 161128_Praktikum_RAJ11_Test_1_layout.csv
