/**
    Some functions involving strings. Added as 'mandatory' utility and to avoid rewritting again and again.

    This is an adapted version of part the code used for analysis in Romero & Ascasibar 2017.

    Mario Romero Spring 2015
**/
#ifndef CONVERSOR_STRING_CPP
#define CONVERSOR_STRING_CPP

#include <sstream> //Only used here
//using namespace std;

string tostring(int a){
	//Conversión entero->string.
	stringstream conv;
	conv<<a;
	return conv.str();
}
string tostring(r_number a){
	//Conversión real->string.
	stringstream conv;
	conv<<a;
	return conv.str();
}
r_number todouble(string A){
	//Conversión string->real
	r_number n;
	stringstream conv;
	conv<<A;
	conv>>n;
	return n;
}
r_number todouble(char *A[]){
//Conversión string->real
	r_number n;
	stringstream conv;
	conv<<A;
	conv>>n;
	return n;
}

bool identify(string word,string phrase){
	//Busca a ver si existe una palabra en un string
	if(phrase.find(word) != string::npos){ //Esta la palabra en la frase.
		return true;}
	return false;
}

string extract(string phrase, string start, string end){
	//Extrae un string entre start (SIN CONTARLO!) y end
	size_t startpos = phrase.find(start) + start.length(); //el problema esq te dara la posición incluyendo el string "start", por lo que hay que añadir la suma que se indica
	size_t endpos = phrase.find(end);
	//size_t length = abs(endpos-startpos); //El valor absoluto esta en el caso de que hayamos puesto las cosas al revés
	size_t length = (endpos-startpos);

	return phrase.substr(startpos,length);

}

#endif
