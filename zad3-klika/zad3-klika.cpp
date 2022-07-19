#include <algorithm>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
using std::cout, std::endl;

struct Sekwencja {
	std::string id;
	std::vector<char> nukleotydy;
	std::vector<int> wiarygodnosc;
	std::vector<int> pozycja_sek_wejsc;
};

std::vector<Sekwencja> wczytywanie_sekwencji(std::ifstream& inputFasta, std::ifstream& inputQual) {
	std::vector<Sekwencja> sekwencje;
	std::string lineFasta;
	std::string lineQual;
	while (std::getline(inputFasta, lineFasta, '>') && std::getline(inputQual, lineQual, '>')) {
		if (lineFasta.empty() || lineQual.empty())
			continue;

		std::string id = lineFasta.substr(0, lineFasta.find(' '));

		lineFasta.erase(0, lineFasta.find('\n'));
		lineFasta.erase(std::remove(lineFasta.begin(), lineFasta.end(), '\n'), lineFasta.end());
		std::vector<char> nukleotydy{ lineFasta.begin(), lineFasta.end() };

		lineQual.erase(0, lineQual.find('\n'));
		std::istringstream numberInteger(lineQual);
		std::vector<int> wiarygodnosc = { std::istream_iterator<int>(numberInteger), std::istream_iterator<int>() };
		std::vector<int> pozycja_sek_wejsc(nukleotydy.size());
		std::iota(std::begin(pozycja_sek_wejsc), std::end(pozycja_sek_wejsc), 1);

		sekwencje.emplace_back(id, nukleotydy, wiarygodnosc, pozycja_sek_wejsc);
	}
	return sekwencje;
}


std::vector<Sekwencja> usuwanie_nukleotydow(int prog_wiarygodnosci, std::vector<Sekwencja> sekwencje) {
	if (prog_wiarygodnosci > 0) {
		for (auto& sekwencja : sekwencje) {
			//cout << sekwencja.id << endl;
			for (int i = 0, j = sekwencja.wiarygodnosc.size(); i < j; ++i) {
				if (sekwencja.wiarygodnosc[i] < prog_wiarygodnosci) {					
					sekwencja.nukleotydy.erase(sekwencja.nukleotydy.begin() + i);
					sekwencja.pozycja_sek_wejsc.erase(sekwencja.pozycja_sek_wejsc.begin() + i);
					sekwencja.wiarygodnosc.erase(sekwencja.wiarygodnosc.begin() + i);
					--i;
					--j;
				}				
			}
		}
		return sekwencje;
	}
	return sekwencje;
}

struct Wierzcholek {
	Sekwencja& sekwencja_wejsciowa;
	size_t pozycja;
	std::string podsekwencja;
};


struct Graf {
	Graf(std::vector<Sekwencja>& sekwencje, size_t dlugosc) {
		for (auto& sekwencja : sekwencje) {

			tworzenie_wierzcholkow(sekwencja, dlugosc);
		}
		//cout << "ok" << endl;
		tworzenie_krawedzi();
	}

	void tworzenie_wierzcholkow(Sekwencja& sekwencja, size_t dlugosc) {

		if (sekwencja.nukleotydy.size() >= dlugosc) {
			//std::cout <<dlugosc << "-literowe podciagi w sekwencji " << sekwencja.id << ":" << std::endl;
			//cout << "Podciag" << "   Pozycja" << endl;
			for (size_t i = 0; i <= sekwencja.nukleotydy.size() - dlugosc; ++i) {
				std::string podsekwencja{};
				for (size_t j = 0; j < dlugosc; ++j) {					
					podsekwencja += sekwencja.nukleotydy[i + j];
				}				
				wierzcholki.emplace_back(sekwencja, sekwencja.pozycja_sek_wejsc[i], podsekwencja);
			}
		}
		else {
			cout << "Brak " << dlugosc << "-literowych podciagow w sekwencji " << sekwencja.id << endl;
		}
	}


	void tworzenie_krawedzi() {
		for (size_t i = 0; i < wierzcholki.size(); ++i) {
			for (size_t j = i + 1; j < wierzcholki.size(); ++j) { 
				if (i != j) {
					if (wierzcholki[i].podsekwencja == wierzcholki[j].podsekwencja) {
						if (&wierzcholki[i].sekwencja_wejsciowa != &wierzcholki[j].sekwencja_wejsciowa) {
							int i_minus_j = wierzcholki[i].pozycja - wierzcholki[j].pozycja;							
							int dziesieciokrotnosc = 10 * wierzcholki[i].podsekwencja.size();
							if (std::abs(i_minus_j) <= dziesieciokrotnosc) {
								std::vector<int> krawedz;
								krawedz.emplace_back(i);
								krawedz.emplace_back(j);
								lista_krawedzi.emplace_back(krawedz);
								krawedz.clear();
							}							
						}
					}
				}
			}
		}
		std::cout << " ilosc krawedzi " << lista_krawedzi.size()<< endl;
	}

	int stopien(int index) {     
		int deg = 0;
		for (int i = 0; i < lista_krawedzi.size(); i++) {
			if (lista_krawedzi[i][0] == index || lista_krawedzi[i][1] == index)
				deg++;      
		}
		return deg;
	}

	std::vector<int> sasiedzi(int index) {     
		std::vector<int> n;
		for (int i = 0; i < lista_krawedzi.size(); i++) {
			if (lista_krawedzi[i][0] == index)
				n.push_back(lista_krawedzi[i][1]);
			if (lista_krawedzi[i][1] == index)
				n.push_back(lista_krawedzi[i][0]);
		}
		return n;
	}

	std::vector<std::vector<int>>  szukanie_kliki(int rozmiar) {
		//cout << "Szukanie kliki o rozmiarze: " << rozmiar << endl;
		if (rozmiar < 2) {
			cout << "Nie znaleziono żadnej kliki!" << endl;
			return{};
		}
		
		std::vector<std::vector<int>> cliques;
		std::vector<int> potential_clique, temp1, temp2;

		for (int i = 0; i < lista_krawedzi.size(); ++i) {
			for (int x = 0; x < 2; ++x) {
				//cout << "deg: " << degree(lista_krawedzi[i][x]) <<endl;
				if (stopien(lista_krawedzi[i][x]) >= rozmiar - 1) {

					potential_clique.push_back(lista_krawedzi[i][x]);    
					temp1 = sasiedzi(lista_krawedzi[i][x]);
					potential_clique.insert(potential_clique.end(), temp1.begin(), temp1.end());
					temp1.clear();
					sort(potential_clique.begin(), potential_clique.end());

					for (int j = 0; j < potential_clique.size(); j++) {      
						temp1.push_back(potential_clique[j]);
						temp2 = sasiedzi(potential_clique[j]);
						temp1.insert(temp1.end(), temp2.begin(), temp2.end());
						temp2.clear();
						sort(temp1.begin(), temp1.end());
						std::vector<int>::iterator it;
						for (int k = 0; k < potential_clique.size(); k++) {				
							if (find(temp1.begin(), temp1.end(), potential_clique[k]) == temp1.end()) {    
								//cout << "usuwanie z potential clique"<<endl;
								potential_clique.erase(find(potential_clique.begin(), potential_clique.end(), potential_clique[k]));
								--k;
							}
						}
						temp1.clear();
					}
					if (potential_clique.size() > rozmiar - 1) {   
						bool already_in = false;
						for (int n = 0; n < cliques.size(); n++) {    
							if (equal(cliques[n].begin(), cliques[n].end(), potential_clique.begin()))
								already_in = true;
						}
						if (!already_in)         
							cliques.push_back(potential_clique);    
					}
					potential_clique.clear();
				}
			}
		}
		if (cliques.empty()) {
			cout << "Nie znaleziono kliki o rozmiarze: " << rozmiar << endl;
			return szukanie_kliki(rozmiar - 1);
		}
		else {
			cout << "Ilosc znalezionych klik o rozmiarze " << rozmiar << ": " << cliques.size() << endl;
			return cliques;
		}
		
	}		
				
	std::vector<Wierzcholek> wierzcholki{};
	std::vector<std::vector<int>> lista_krawedzi{};
};


int main() {
	std::string plik;
	int prog_wiarygodnosci, dlugosc;

	std::cout << "Podaj nazwe pliku: ";
	std::cin >> plik;

	std::cout << "Podaj dlugosc podciagu: (miedzy 4 a 9) ";
	std::cin >> dlugosc;

	std::cout << "Podaj prog wiarygodnosci: ";
	std::cin >> prog_wiarygodnosci;

	std::ifstream fasta("C:\\Users\\Antek\\Documents\\akwb-lab\\zad3\\fasta\\" + plik + ".fasta");
	std::ifstream qual("C:\\Users\\Antek\\Documents\\akwb-lab\\zad3\\qual\\" + plik + ".qual");

	auto sekwencje = wczytywanie_sekwencji(fasta, qual);
	
	auto nowe_sekwencje = usuwanie_nukleotydow(prog_wiarygodnosci, sekwencje);
	
	//cout << nowe_sekwencje.size();

	Graf graf(nowe_sekwencje, dlugosc);
	std::cout << "Instancja: " << plik << "\n";
	std::cout << "Dlugosc: " << dlugosc << "\t Prog wiarygodnosci: " << prog_wiarygodnosci << "\n";

	int index;
	std::vector<std::vector<int>> kliki = graf.szukanie_kliki(sekwencje.size());
	if (kliki.empty()) {
		cout << "Nie znaleziono motywu o podanej dlugosci we wszystkich sekwencjach instancji" << endl;
	}
	else {		
		for (int i = 0; i < kliki.size(); i++) {
			index = kliki[i][0];
			cout << endl;			
			cout << "Podciag :  " << graf.wierzcholki[index].podsekwencja << "\n";
			for (int j = 0; j < kliki[i].size(); j++) {
				index = kliki[i][j];
				cout << "\t " << graf.wierzcholki[index].sekwencja_wejsciowa.id << "  Pozycja: " << graf.wierzcholki[index].pozycja << "\n";
			}
		}
	}

	return 0;
}