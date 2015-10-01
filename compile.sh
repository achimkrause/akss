gcc smith.c matrix.c abelian.c test.c main.c -std=c99 -lgmp -o main 2> log

if [[ -s log ]]; then 
less log
fi;
