# ZPSplanetdet
Repozytorium Zespołowego Projektu Studenckiego, realizowanego w ramach studiów magisterskich na Wydziale Fizyki UW. 

# TO DO:

| Status | Task | Person | Notes |
|--------|------|--------|-------|
|* | find out how degenerate cases were treated in Suzuki paper|Gabriela| |
|* | calculate N\_exp in simplest case analytically |Gabriela| |
| | MAYBE: speed-up N\_exp calculations | | this task is useful, but not required |
|* | write code that repeats the analysis of Suzuki+16 (everything except N\_exp calculation) | Pawel, Maciej | prepared likelihood and probality functions |
\ | extend Suzuki+16 analysis to better include degenerate solutions | | |
\ | extend Suzuki+16 analysis to include uncertainties of s and q for each planet | | | 

Status (- nothing done, * in progress, x done)

# ALREADY DONE:

| Status | Task | Person | Notes |
|--------|------|--------|-------|
|x| add papers used in Suzuki| Pawel | MOA-2012-BLG-355 paper missing |
|x| calculate f(s1,q)/f(s2,q) for degenerate cases|Maciej| |
|x| find and download raw S(s, q) data from Suzuki+16|RP| |
|x| plot S(s, q), i.e., fig 6 from paper|Marcin| |
|x| add raw planet data (s,q) | Pawel | errors should be added next |
|x| implement function calculating N\_exp|Pawel, Marcin, Jakub| Pawel's solution implemented with target *f* function |
|x| test N\_exp calculation using S(s, q) and f(s, q) from S+16| | |



# zasady pisania kodu
- nazwy zmiennych i funkcji są znaczące i określają po co jest dana zmienna/funkcja
- nazwy funkcji zaczynają się od czasownika
- pod komendą def jest docstring, który opisuje funkcję i ewentualnie wejście oraz wyjście
- wcięcia robimy czterema spacjami, a nie tabulatorem
