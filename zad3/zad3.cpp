#include <algorithm>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Sekwencja {
	std::string id;
	std::vector<char> nukleotydy;
	std::vector<int> wiarygodnosc;
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

		sekwencje.emplace_back(id, nukleotydy, wiarygodnosc);
	}
	return sekwencje;
}


struct Wierzcholek {
	Sekwencja& sekwencja_wejsciowa;
	size_t pozycja;
	std::string podsekwencja;
};


struct Graf {
	Graf(std::vector<Sekwencja>& sekwencja, size_t dlugosc, size_t prog_wiarygodnosci) {
		for (auto& each : sekwencja)
			tworzenie_wierzcholkow(each, dlugosc, prog_wiarygodnosci);
		tworzenie_krawedzi();
	}


	void tworzenie_wierzcholkow(Sekwencja& sekwencja, size_t dlugosc, size_t prog_wiarygodnosci) {
		for (size_t i = 0; i < sekwencja.nukleotydy.size(); ++i) {
			if (sekwencja.wiarygodnosc[i] < prog_wiarygodnosci)
				continue;

			std::string podsekwencja{};
			for (size_t j = 0, left = dlugosc; j < left; ++j) {
				if (sekwencja.nukleotydy.size() - i - j < left - j)
					return;

				if (sekwencja.wiarygodnosc[i + j] < prog_wiarygodnosci)
					++left;
				else

					podsekwencja += sekwencja.nukleotydy[i + j];
			}
			wierzcholki.emplace_back(sekwencja, i, podsekwencja);
		}
	}


	void tworzenie_krawedzi() {
		for (size_t i = 0; i < wierzcholki.size(); ++i)
			for (size_t j = 0; j < wierzcholki.size(); ++j)
				if (i != j)
					if (wierzcholki[i].podsekwencja == wierzcholki[j].podsekwencja)
						if (&wierzcholki[i].sekwencja_wejsciowa != &wierzcholki[j].sekwencja_wejsciowa)
							if (wierzcholki[i].pozycja - wierzcholki[j].pozycja < 10 * wierzcholki[i].podsekwencja.size() || wierzcholki[j].pozycja - wierzcholki[i].pozycja < 10 * wierzcholki[j].podsekwencja.size())//czy różnica w pozycji nie jest większa od 10krotnej długości podciągu
								krawedzie[&wierzcholki[i]].push_back(&wierzcholki[j]);

	}


	std::vector<Wierzcholek*> cliqueGenerator() {
		Wierzcholek* biezacy_wierzcholek = nullptr;
		for (const auto& [wierzcholek, wierzcholki_sasiadujace] : krawedzie)
			if (wierzcholki_sasiadujace.size() > krawedzie[biezacy_wierzcholek].size())
				biezacy_wierzcholek = wierzcholek;

		if (!biezacy_wierzcholek)
			return {};

		std::vector<Sekwencja*> sequences{ &biezacy_wierzcholek->sekwencja_wejsciowa };
		std::vector<Wierzcholek*> clique{ biezacy_wierzcholek };

		bool finished = false;
		while (!finished) {
			finished = true;

			for (const auto& wierzcholek_sasiadujacy : krawedzie[clique.back()]) {
				if (std::find(clique.begin(), clique.end(), wierzcholek_sasiadujacy) != clique.end())
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


				if (std::find(sequences.begin(), sequences.end(), &wierzcholek_sasiadujacy->sekwencja_wejsciowa) != sequences.end())
					continue;

				clique.push_back(wierzcholek_sasiadujacy);
				sequences.push_back(&wierzcholek_sasiadujacy->sekwencja_wejsciowa);
				finished = false;
				break;
			}
		}
		return clique;
	}

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
	Graf graf(sekwencje, dlugosc, prog_wiarygodnosci);

	std::cout << "plik: " << plik << " dlugosc: " << dlugosc << " prog wiarygodnosci: " << prog_wiarygodnosci << "\n";
	for (const auto& vertex : graf.cliqueGenerator())
		std::cout << "ID sekwencji: " << vertex->sekwencja_wejsciowa.id
		<< " Pozycja: " << vertex->pozycja
		<< " Nukleotyd: " << vertex->podsekwencja << "\n";

	return 0;
}