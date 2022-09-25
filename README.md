# ZPSplanetdet
Repozytorium Zespołowego Projektu Studenckiego, realizowanego w ramach studiów magisterskich na Wydziale Fizyki UW. 
Plik raportu na OverLeaf: https://www.overleaf.com/read/gbpfvgwbykdc

# TO DO:

| Status | Task | Person | Notes |
|--------|------|--------|-------|
|\* | extend Suzuki+16 analysis to better include degenerate solutions | Maciej, Gabriela | |
| | prepare text of the report | | |
| | prepare plot for the report | | |
| | MAYBE: speed-up N\_exp calculations | | this task is useful, but not required |
| | MAYBE: combine common (the same or almost the same) functions from mcmc\_kuba.py and mcmc\_degen.py | | |
| | MAYBE: print results to screen (similar to corner call)
| | MAYBE: extend Suzuki+16 analysis to include uncertainties of s and q for each planet | | | 

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
|x| find out how degenerate cases were treated in Suzuki paper|Gabriela| [degenerate\_cases.txt](degenerate_cases.txt) |
|x| interpolation of S(s,q)|Jakub| |
|\*| calculate N\_exp in simplest case analytically |Gabriela| no need to finish it|
|x| calculate weights based on Dchi^2|Gabriela, Maciej| |
|x| write code that repeats the analysis of Suzuki+16 | Pawel, Maciej, Jakub, Gabriela, Marcin | |
|x| basic output of results | Jakub | |
|x| write code that repeats the analysis of Suzuki+16 | Pawel, Maciej, Jakub, Gabriela, Marcin | |
|x| speed-up code by reducing interpolation function calls | Maciej, Gabriela | |




# zasady pisania kodu
- nazwy zmiennych i funkcji są znaczące i określają po co jest dana zmienna/funkcja
- nazwy funkcji zaczynają się od czasownika
- pod komendą def jest docstring, który opisuje funkcję i ewentualnie wejście oraz wyjście
- wcięcia robimy czterema spacjami, a nie tabulatorem
