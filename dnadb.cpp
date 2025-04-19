// CMSC 341 - Spring 25 - Project 4
#include "dnadb.h"
DnaDb::DnaDb(int size, hash_fn hash, prob_t probing = DEFPOLCY){
    
}

DnaDb::~DnaDb(){
    
}

void DnaDb::changeProbPolicy(prob_t policy){
    
}

bool DnaDb::insert(DNA dna){
    
}

bool DnaDb::remove(DNA dna){
    
}

const DNA DnaDb::getDNA(string sequence, int location) const{
    
}

bool DnaDb::updateLocId(DNA dna, int location){
    
}

float DnaDb::lambda() const {
      
}

float DnaDb::deletedRatio() const {
    
}

void DnaDb::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}