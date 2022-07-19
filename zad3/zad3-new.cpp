#include <algorithm>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Sequence {
	std::string sequenceID;
	std::vector<char> nucleotides;
	std::vector<int> veracity;
};

std::vector<Sequence> sequenceParser(std::ifstream& inputFasta, std::ifstream& inputQual) {
	std::vector<Sequence> sequences;
	std::string lineFasta;
	std::string lineQual;
	while (std::getline(inputFasta, lineFasta, '>') && std::getline(inputQual, lineQual, '>')) {
		if (lineFasta.empty() || lineQual.empty())
			continue;

		std::string sequenceID = lineFasta.substr(0, lineFasta.find(' '));

		lineFasta.erase(0, lineFasta.find('\n'));
		lineFasta.erase(std::remove(lineFasta.begin(), lineFasta.end(), '\n'), lineFasta.end());
		std::vector<char> nucleotides{ lineFasta.begin(), lineFasta.end() };

		lineQual.erase(0, lineQual.find('\n'));
		std::istringstream numberInteger(lineQual);
		std::vector<int> veracity = { std::istream_iterator<int>(numberInteger), std::istream_iterator<int>() };

		sequences.emplace_back(sequenceID, nucleotides, veracity);
	}
	return sequences;
}


struct Peak {
	Sequence& sequence;
	size_t position;
	std::string nucleotides;
};

struct Graph {
	Graph(std::vector<Sequence>& sequences, size_t length, size_t threshold) {
		for (auto& each : sequences)
			peakGenerator(each, length, threshold);
		edgesGenerator();
	}

	void peakGenerator(Sequence& sequence, size_t length, size_t threshold) {
		for (size_t i = 0; i < sequence.nucleotides.size(); ++i) {
			if (sequence.veracity[i] < threshold)
				continue;

			std::string subSequence{};
			for (size_t j = 0, left = length; j < left; ++j) {
				if (sequence.nucleotides.size() - i - j < left - j)
					return;

				if (sequence.veracity[i + j] < threshold)
					++left;
				else

					subSequence += sequence.nucleotides[i + j];
			}
			peaks.emplace_back(sequence, i, subSequence);
		}
	}


	// void edgesGenerator() {
	// 	for (size_t i = 0; i < peaks.size(); ++i)
	// 		for (size_t j = 0; j < peaks.size(); ++j)
	// 			if (i != j)
	// 				if (peaks[i].nucleotides == peaks[j].nucleotides)
	// 					if (&peaks[i].sequence != &peaks[j].sequence)
	// 						edges[&peaks[i]].push_back(&peaks[j]);
	// }


	void edgesGenerator() {
		for (size_t i = 0; i < peaks.size(); ++i)
			for (size_t j = 0; j < peaks.size(); ++j)
				if (i != j)
					if (peaks[i].nucleotides == peaks[j].nucleotides)
						if (&peaks[i].sequence != &peaks[j].sequence)
							if (peaks[i].position - peaks[j].position < 10* peaks[i].nucleotides.size() || peaks[j].position-peaks[i].position < 10*peaks[j].nucleotides.size())//czy różnica w pozycji nie jest większa od 10krotnej długości podciągu
								edges[&peaks[i]].push_back(&peaks[j]);
							
	}

	std::vector<Peak*> cliqueGenerator() {
		Peak* peakCurrent = nullptr;
		for (const auto& [peak, peakNext] : edges)
			if (peakNext.size() > edges[peakCurrent].size())
				peakCurrent = peak;


		if (!peakCurrent)
			return {};

		std::vector<Sequence*> sequences{ &peakCurrent->sequence };
		std::vector<Peak*> clique{ peakCurrent };

		bool finished = false;
		while (!finished) {
			finished = true;

			for (const auto& peakNext : edges[clique.back()]) {
				if (std::find(clique.begin(), clique.end(), peakNext) != clique.end())
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


				if (std::find(sequences.begin(), sequences.end(), &peakNext->sequence) != sequences.end())
					continue;

				clique.push_back(peakNext);
				sequences.push_back(&peakNext->sequence);
				finished = false;
				break;
			}
		}
		return clique;
	}


	std::vector<Peak> peaks{};
	std::map<Peak*, std::vector<Peak*>> edges{};
};


int main() {

	std::string filename;
	int length;
	int threshold;

	std::cout << "Podaj nazwe pliku: ";
	std::cin >> filename;

	std::cout << "Podaj dlugosc podciagu: (miedzy 4 a 9) ";
	std::cin >> length;

	std::cout << "Podaj prog wiarygodnosci: ";
	std::cin >> threshold;

	std::ifstream inputFasta("C:\\Users\\Antek\\Documents\\akwb-lab\\projekt3\\fasta\\" + filename + ".fasta");
	std::ifstream inputQual("C:\\Users\\Antek\\Documents\\akwb-lab\\projekt3\\qual\\" + filename + ".qual");
	std::cout << "C:\\Users\\Antek\\Documents\\akwb - lab\\projekt3\\qual\\" + filename + ".qual";

	auto sequences = sequenceParser(inputFasta, inputQual);
	Graph graph(sequences, length, threshold);

	std::cout << "file: " << filename << " length: " << length << " threshold: " << threshold << "\n";
	for (const auto& vertex : graph.cliqueGenerator())
		std::cout << "SequenceID: " << vertex->sequence.sequenceID
		<< " Position: " << vertex->position
		<< " Nucleotide: " << vertex->nucleotides << "\n";

	return 0;
}