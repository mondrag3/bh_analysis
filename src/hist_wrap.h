#ifndef hist_wrap_h
#define hist_wrap_h

#include <iostream>
#include <string>
#include <vector>
#include <map>

class TH1F;

class hist {
  struct binning {
    int nbins; double min, max;
    binning();
    ~binning();
  };
  static std::map<std::string,binning> binnings;
  static std::vector<const hist*> all;
  static std::string binning_name_regex_pattern;
  static binning default_binning;

  const binning get_binning(const std::string& hist_name);

  const binning b;
  std::pair<int,double> underflow, overflow;
  TH1F * h;

public:
  hist();
  hist(const std::string& name,const std::string& title);
  ~hist();

  void Fill(double x, double w=1.);
  void FillOverflow(double weight);
  static void read_binnings(const char* filename, const char* regex="");
  static std::ostream& print_overflow(std::ostream& out=std::cout);
  static void delete_all();
};

#endif
