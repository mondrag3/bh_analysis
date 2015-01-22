#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "rapidxml-1.13/rapidxml.hpp"

using namespace std;
using namespace rapidxml;

struct particle {
  string name;
  int pid;
  vector<particle> daughters;
  particle(const xml_node<> *node)
  : name(node->first_attribute("name")->value()),
    pid([](const xml_attribute<char> *attr) {
      return (attr ? atoi(attr->value()) : 0);
    }(node->first_attribute("pid")) )
  {
    for (xml_node<> *p = node->first_node(); p; p = p->next_sibling())
      daughters.emplace_back(p);
  }
};

int main(int argc, char **argv)
{
  if (argc!=2) {
    cout << "Usage: " << argv[0] << " file.xml" << endl;
    exit(0);
  }

  cout << "XML file: " << argv[1] << endl;

  xml_document<> doc;
	// Read the xml file into a vector
	ifstream theFile(argv[1]);
	vector<char> buffer((istreambuf_iterator<char>(theFile)),
                      istreambuf_iterator<char>());
	buffer.push_back('\0');
	// Parse the buffer using the xml file parsing library into doc
	doc.parse<0>(buffer.data());
  // Find root node
	const xml_node<> *analysis   = doc.first_node("analysis");
  const xml_node<> *particles  = analysis->first_node("particles");
  // const xml_node<> *jets       = analysis->first_node("jets");
  // const xml_node<> *histograms = analysis->first_node("histograms");

  for (xml_node<> *_p = particles->first_node(); _p; _p = _p->next_sibling()) {
    particle p(_p);
    cout << p.name << ' ' << p.pid << endl;
    for (auto& pp : p.daughters)
      cout << "  " << pp.name << ' ' << pp.pid << endl;
	}

  return 0;
}
