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

//std::vector<Sekwencja> usuwanie_nukleotydow(int prog_wiarygodnosci, std::vector<Sekwencja> sekwencje) {
//	if (prog_wiarygodnosci > 0) {
//		std::vector<Sekwencja> nowe_sekwencje;
//		for (auto& sekwencja : sekwencje) {
//			std::vector<char> temp_nukleotydy;
//			std::vector<int> temp_wiarygodnosc;
//			std::vector<int> temp_pozycja;
//			for (size_t i = 0, i < sekwencja.wiarygodnosc.size(); ++i) {
//				if (sekwencja.wiarygodnosc[i] >= prog_wiarygodnosci) {
//					temp_nukleotydy.emplace_back(sekwencja.nukleotydy[i]);
//					temp_wiarygodnosc.emplace_back(sekwencja.wiarygodnosc[i]);
//					temp_pozycja.emplace_back(sekwencja.pozycja_sek_wejsc[i]);
//				}
//			}
//			nowe_sekwencje.emplace_back(sekwencja.id, temp_nukleotydy, temp_wiarygodnosc, temp_pozycja);
//		}
//		return nowe_sekwencje;
//	}
//	else {
//		return sekwencje;
//	}
//}
std::vector<Sekwencja> usuwanie_nukleotydow(int prog_wiarygodnosci, std::vector<Sekwencja> sekwencje) {
	if (prog_wiarygodnosci > 0) {		
		for (auto& sekwencja : sekwencje) {	
			//cout << sekwencja.id << endl;
			for (int i = 0, j = sekwencja.wiarygodnosc.size(); i < j; ++i) {
				if (sekwencja.wiarygodnosc[i] < prog_wiarygodnosci) {
					//std::cout << "mniejsza: "<< sekwencja.wiarygodnosc[i] << endl;
					sekwencja.nukleotydy.erase(sekwencja.nukleotydy.begin() + i);
					sekwencja.pozycja_sek_wejsc.erase(sekwencja.pozycja_sek_wejsc.begin()+i);
					sekwencja.wiarygodnosc.erase(sekwencja.wiarygodnosc.begin() + i);
					i--;
					j--;
				}
				else {
					//std::cout << sekwencja.wiarygodnosc[i] <<" ";					
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
		cout << "ok"<<endl;			
		tworzenie_krawedzi();
	}

	void tworzenie_wierzcholkow(Sekwencja& sekwencja, size_t dlugosc) {
		
		if (sekwencja.nukleotydy.size() >= dlugosc) {
			//std::cout <<dlugosc << "-literowe podciagi w sekwencji " << sekwencja.id << ":" << std::endl;
			//cout << "Podciag" << "   Pozycja" << endl;
			for (size_t i = 0; i <= sekwencja.nukleotydy.size() - dlugosc; ++i) {			
				std::string podsekwencja{};
				for (size_t j = 0; j < dlugosc; ++j) {		
					/*cout << "i: " << i;
					cout << " j: " << j;
					cout << " ";*/
					podsekwencja += sekwencja.nukleotydy[i + j];				
				}

				//cout << podsekwencja << "      " ;
				//cout << sekwencja.pozycja_sek_wejsc[i] << endl;
				wierzcholki.emplace_back(sekwencja, sekwencja.pozycja_sek_wejsc[i], podsekwencja);
			}
		}
		else {
			cout << "Brak " << dlugosc << "-literowych podciagow w sekwencji " << sekwencja.id << endl;
		}
	}	


	void tworzenie_krawedzi() {
		for (size_t i = 0; i < wierzcholki.size(); ++i){
			for (size_t j = i+1; j < wierzcholki.size(); ++j){ //for (size_t j = 0; j < wierzcholki.size(); ++j){
				if (i != j) {				
					if (wierzcholki[i].podsekwencja == wierzcholki[j].podsekwencja){
						if (&wierzcholki[i].sekwencja_wejsciowa != &wierzcholki[j].sekwencja_wejsciowa) {
							int i_minus_j = wierzcholki[i].pozycja - wierzcholki[j].pozycja;
							/*int j_minus_i = wierzcholki[j].pozycja - wierzcholki[i].pozycja;*/
							int dziesieciokrotnosc = 10 * wierzcholki[i].podsekwencja.size();							
							if (std::abs(i_minus_j) <= dziesieciokrotnosc /*&& std::abs(j_minus_i) <= dziesieciokrotnosc*/) {
								krawedzie[&wierzcholki[i]].push_back(&wierzcholki[j]);
								//std::cout << "dziesieciokrotnosc: " << dziesieciokrotnosc;
								//std::cout << " roznica: " << std::abs(i_minus_j) << endl;
							}								
						}	
					}
				}
			}
		}
	}


	std::vector<std::vector<Wierzcholek*>>  szukanie_kliki(int rozmiar) {
		cout << "Szukanie kliki o rozmiarze: " << rozmiar << endl;
		if (rozmiar < 2) {
			cout << "Nie znaleziono żadnej kliki!"  << endl;
			return{};
		}
		std::vector<Wierzcholek*> wierzcholki_max_deg;
		Wierzcholek* biezacy_wierzcholek = nullptr;

		for (const auto& [wierzcholek, wierzcholki_sasiadujace] : krawedzie)
			if (wierzcholki_sasiadujace.size() >= rozmiar -1)
				wierzcholki_max_deg.emplace_back(wierzcholek);

		if (wierzcholki_max_deg.empty()){
			cout << " (max-deg) Nie znaleziono kliki o rozmiarze " << rozmiar << endl;
			return szukanie_kliki(rozmiar-1);
			return {};
		}
		std::vector<std::vector<Wierzcholek*>> kliki;
		for (const auto wierzcholek_max_deg : wierzcholki_max_deg) {
			biezacy_wierzcholek = wierzcholek_max_deg;
			std::vector<Sekwencja*> sekwencje{ &biezacy_wierzcholek->sekwencja_wejsciowa };
			std::vector<Wierzcholek*> klika{ biezacy_wierzcholek };
			bool finished = false;
			while (!finished) {
				finished = true;

				for (const auto& wierzcholek_sasiadujacy : krawedzie[klika.back()]) {
					if (std::find(klika.begin(), klika.end(), wierzcholek_sasiadujacy) != klika.end())
						continue;
					/*
					int counter = 1;
					for (const auto& peakFirst : edges[clique.back()])
					{
						if (peakFirst != clique.end())
							continue;
						for (const auto& peakNext : edges[sequences.end()])
						{

							if (peakNext == sequences.end())
								continue;
							if (clique.begin() == peakNext)
							{
								counter += 1;
							}
							if (counter == peakNext)
							{
								break;
							}
					*/
					if (std::find(sekwencje.begin(), sekwencje.end(), &wierzcholek_sasiadujacy->sekwencja_wejsciowa) != sekwencje.end())
						continue;

					klika.push_back(wierzcholek_sasiadujacy);
					sekwencje.push_back(&wierzcholek_sasiadujacy->sekwencja_wejsciowa);
					finished = false;
					break;
				}
			}
			if (klika.size() == rozmiar ) {	
				//if (kliki.empty()) {
				//	std::sort(klika.begin(), klika.end());
				//	kliki.emplace_back(klika);
				//}
				//else if (std::find(kliki.begin(), kliki.end(), klika) == kliki.end()) {
				//	//if (  klika not in kliki  )
				//	std::sort(klika.begin(), klika.end());
				//	kliki.emplace_back(klika);
				//}			
				//std::sort(klika.begin(), klika.end());
				kliki.emplace_back(klika);
			}	
		}
		if (kliki.empty()) {
			cout << "Nie znaleziono kliki o rozmiarze: " << rozmiar << endl;
			return szukanie_kliki(rozmiar - 1);
		}		
		/*cout <<"rozmiar przed unique: " << kliki.size()<<endl;
		std::vector<std::vector<Wierzcholek*>>::iterator ip;
		std::sort(kliki.begin(), kliki.end());
		ip = std::unique(kliki.begin(), kliki.end());
		kliki.resize(std::distance(kliki.begin(), ip));
		cout<<"rozmiar po unique: " << kliki.size()<<endl;*/
		cout << "Ilosc znalezionych klik o rozmiarze " << rozmiar << ": " << kliki.size() << endl;

		return kliki;
	}


	//std::vector<Wierzcholek*> cliqueGenerator() {
	//	Wierzcholek* biezacy_wierzcholek = nullptr;
	//	for (const auto& [wierzcholek, wierzcholki_sasiadujace] : krawedzie)
	//		if (wierzcholki_sasiadujace.size() > krawedzie[biezacy_wierzcholek].size())
	//			biezacy_wierzcholek = wierzcholek;
	//	if (!biezacy_wierzcholek)
	//		return {};
	//	std::vector<Sekwencja*> sekwencje{ &biezacy_wierzcholek->sekwencja_wejsciowa };
	//	std::vector<Wierzcholek*> klika{ biezacy_wierzcholek };
	//	bool finished = false;
	//	while (!finished) {
	//		finished = true;
	//		for (const auto& wierzcholek_sasiadujacy : krawedzie[klika.back()]) {
	//			if (std::find(klika.begin(), klika.end(), wierzcholek_sasiadujacy) != klika.end())
	//				continue;
	//			/*
	//			int counter = 1;
	//			for (const auto& peakFirst : edges[clique.back()])
	//			{
	//				if (peakFirst != clique.end())
	//					continue;
	//				for (const auto& peakNext : edges[sequences.end()])
	//				{
	//					if (peakNext == sequences.end())
	//						continue;
	//					if (clique.begin() == peakNext)
	//					{
	//						counter += 1;
	//					}
	//					if (counter == peakNext)
	//					{
	//						break;
	//					}
	//			*/
	//			if (std::find(sekwencje.begin(), sekwencje.end(), &wierzcholek_sasiadujacy->sekwencja_wejsciowa) != sekwencje.end())
	//				continue;
	//			klika.push_back(wierzcholek_sasiadujacy);
	//			sekwencje.push_back(&wierzcholek_sasiadujacy->sekwencja_wejsciowa);
	//			finished = false;
	//			break;
	//		}
	//	}
	//	return klika;
	//}

	std::vector<Wierzcholek> wierzcholki{};
	std::map<Wierzcholek*, std::vector<Wierzcholek*>> krawedzie{};
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

	

	//for (auto& sekwencja : sekwencje) {
	//	std::cout << sekwencja.id << std::endl;
	//	for (auto& pozycja : sekwencja.pozycja_sek_wejsc)
	//		std::cout << pozycja << " ";
	//	std::cout << std::endl;
	//	for (auto& nukleotyd : sekwencja.nukleotydy)
	//		std::cout << nukleotyd << " ";
	//	std::cout << std::endl;
	//	/*for (auto& wiarygodnosc : sekwencja.wiarygodnosc)
	//		std::cout << wiarygodnosc << " ";*/
	//	std::cout << std::endl;
	//}


	auto nowe_sekwencje = usuwanie_nukleotydow(prog_wiarygodnosci, sekwencje);

	/*std::cout << "PO ZMIANIE: " << std::endl;
	std::cout << std::endl;*/
	//for (auto& sekwencja : nowe_sekwencje) {
	//	std::cout << sekwencja.id << std::endl;
	//	for (auto& pozycja : sekwencja.pozycja_sek_wejsc)
	//		std::cout << pozycja << " ";
	//	std::cout << std::endl;
	//	for (auto& nukleotyd : sekwencja.nukleotydy)
	//		std::cout << nukleotyd << " ";
	//	std::cout << std::endl;
	//	/*for (auto& wiarygodnosc : sekwencja.wiarygodnosc)
	//		std::cout << wiarygodnosc <<" ";*/
	//	std::cout << std::endl;
	//}
	cout << nowe_sekwencje.size();

	Graf graf(nowe_sekwencje, dlugosc);
	
	std::cout << "plik: " << plik << " dlugosc: " << dlugosc << " prog wiarygodnosci: " << prog_wiarygodnosci << "\n";


	std::vector<std::vector<Wierzcholek*>> kliki = graf.szukanie_kliki(sekwencje.size());


	for (const auto& klika : kliki) {
		std::cout << "Klika:  " << klika[0]->podsekwencja << endl;
		for (const auto& vertex : klika) {
			std::cout << "ID sekwencji: " << vertex->sekwencja_wejsciowa.id;
			std::cout << " Pozycja: " << vertex->pozycja << endl;
			//std::cout << " Podciag: " << vertex->podsekwencja << endl;
		}
		std::cout << endl;
	}
	//for (const auto& klika : graf.cliqueGenerator(sekwencje.size())) {
	//	std::cout<< "Klika:  " << klika[0]->podsekwencja << endl;
	//	for (const auto& vertex : klika){
	//		std::cout << "ID sekwencji: " << vertex->sekwencja_wejsciowa.id;
	//		std::cout << " Pozycja: " << vertex->pozycja << endl;
	//		//std::cout << " Podciag: " << vertex->podsekwencja << endl;
	//	}
	//	std::cout << endl;
	//}	
	return 0;
}