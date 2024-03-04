#ifndef B1ScintHit_h
#define B1ScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class B1ScintHit : public G4VHit
{
  public:
    B1ScintHit();
    B1ScintHit(const B1ScintHit&);
    virtual ~B1ScintHit();

    // operators
    const B1ScintHit& operator=(const B1ScintHit&);
    G4bool operator==(const B1ScintHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;
      
  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the  sensitive volume
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using B1ScintHitsCollection = G4THitsCollection<B1ScintHit>;

extern G4ThreadLocal G4Allocator<B1ScintHit>* B1ScintHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* B1ScintHit::operator new(size_t)
{
  if (!B1ScintHitAllocator) {
    B1ScintHitAllocator = new G4Allocator<B1ScintHit>;
  }
  void *hit;
  hit = (void *) B1ScintHitAllocator->MallocSingle();
  return hit;
}

inline void B1ScintHit::operator delete(void *hit)
{
  if (!B1ScintHitAllocator) {
    B1ScintHitAllocator = new G4Allocator<B1ScintHit>;
  }
  B1ScintHitAllocator->FreeSingle((B1ScintHit*) hit);
}

inline void B1ScintHit::Add(G4double de, G4double dl) {
  fEdep += de; 
  fTrackLength += dl;
}

inline G4double B1ScintHit::GetEdep() const { 
  return fEdep; 
}

inline G4double B1ScintHit::GetTrackLength() const { 
  return fTrackLength; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
