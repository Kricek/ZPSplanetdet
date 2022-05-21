# ZPSplanetdet
Repozytorium Zespołowego Projektu Studenckiego, realizowanego w ramach studiów magisterskich na Wydziale Fizyki UW. 

# TO DO

| Status | Task | Person | Notes |
|--------|------|--------|-------|
|x| add papers used in Suzuki| Pawel | MOA-2012-BLG-355 paper missing |
|x| calculate f(s1,q)/f(s2,q) for degenerate cases|Maciej| |
|* | find out how degenerate cases were treated in Suzuki paper|Gabriela| |
|x| find and download raw S(s, q) data from Suzuki+16|RP| |
|x| plot S(s, q), i.e., fig 6 from paper|Marcin| |
|* | calculate N\_exp in simplest case analytically |Gabriela| |
|* | implement function calculating N\_exp|Pawel, Marcin| Pawel's solution implemented with target *f* function |
|-| test N\_exp calculation using S(s, q) and f(s, q) from S+16| | |
|x| add raw planet data (s,q) | Pawel | errors should be added next |
|* | write code that repeats the analysis of Suzuki+16 (everything except N\_exp calculation) | Pawel | prepared likelihood and probality functions |

Status (- nothing done, * in progress, x done)

# zasady pisania kodu
- nazwy zmiennych i funkcji są znaczące i określają po co jest dana zmienna/funkcja
- nazwy funkcji zaczynają się od czasownika
- pod komendą def jest docstring, który opisuje funkcję i ewentualnie wejście oraz wyjście
- wcięcia robimy czterema spacjami, a nie tabulatorem
