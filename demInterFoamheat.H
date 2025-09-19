#include "demInterFoamBaseheat.H"

class demInterFoamheat : public demInterFoamBaseheat {
public:
  demInterFoamheat();
  virtual ~demInterFoamheat();
  virtual void run(double t);
  
  pimpleControl *pimple_;
};

