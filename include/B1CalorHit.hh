#ifndef B1CalorHit_h
#define B1CalorHit_h 1

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

class B1CalorHit : public G4VHit
{
  public:
    B1CalorHit();
    B1CalorHit(const B1CalorHit&);
    virtual ~B1CalorHit();

    // operators
    const B1CalorHit& operator=(const B1CalorHit&);
    G4bool operator==(const B1CalorHit&) const;

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

using B1CalorHitsCollection = G4THitsCollection<B1CalorHit>;

extern G4ThreadLocal G4Allocator<B1CalorHit>* B1CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* B1CalorHit::operator new(size_t)
{
  if (!B1CalorHitAllocator) {
    B1CalorHitAllocator = new G4Allocator<B1CalorHit>;
  }
  void *hit;
  hit = (void *) B1CalorHitAllocator->MallocSingle();
  return hit;
}

inline void B1CalorHit::operator delete(void *hit)
{
  if (!B1CalorHitAllocator) {
    B1CalorHitAllocator = new G4Allocator<B1CalorHit>;
  }
  B1CalorHitAllocator->FreeSingle((B1CalorHit*) hit);
}

inline void B1CalorHit::Add(G4double de, G4double dl) {
  fEdep += de; 
  fTrackLength += dl;
}

inline G4double B1CalorHit::GetEdep() const { 
  return fEdep; 
}

inline G4double B1CalorHit::GetTrackLength() const { 
  return fTrackLength; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
