// CMSC 341 - Spring 25 - Project 4
#ifndef DNADB_H
#define DNADB_H
#include <iostream>
#include <string>
#include "math.h"
using namespace std;
class Grader;   // forward declaration, will be used for grdaing
class Tester;   // forward declaration, will be used for testing
class DNA;      // forward declaration
class DnaDb;    // forward declaration
const int MINPRIME = 101;   // Min size for hash table
const int MAXPRIME = 99991; // Max size for hash table
const int MINLOCID = 100000;// Min Location ID
const int MAXLOCID = 999999;// Max Location ID
typedef unsigned int (*hash_fn)(string);     // declaration of hash function
enum prob_t {QUADRATIC, DOUBLEHASH, LINEAR}; // types of collision handling policy
#define DEFPOLCY QUADRATIC
const int MAX = 4;
const char ALPHA[MAX] = {'A', 'C', 'G', 'T'};
class DNA{
    public:
    friend class Grader;
    friend class Tester;
    friend class DnaDb;
    DNA(string sequence="", int location=0, bool used=false){
        m_sequence=sequence; m_location=location; m_used=used;
    }
    string getSequence() const {return m_sequence;}
    int getLocId() const {return m_location;}
    bool getUsed() const {return m_used;}
    void setSequence(string seq) {m_sequence=seq;}
    void setLocID(int id) {m_location=id;}
    void setUsed(bool used) {m_used=used;}
    // the following function is a friend function
    friend ostream& operator<<(ostream& sout, const DNA *dna ){
        if ((dna != nullptr) && !(dna->getSequence().empty()))
            sout << dna->getSequence() << " (" << dna->getLocId() << ", "<< dna->getUsed() <<  ")";
        else
            sout << "";
        return sout;
    }
    // the following function is a friend function
    friend bool operator==(const DNA& lhs, const DNA& rhs){
        // since the uniqueness of an object is defined by sequence and location ID
        // the equality operator considers only those two criteria
        return ((lhs.getSequence() == rhs.getSequence()) && (lhs.getLocId() == rhs.getLocId()));
    }
    // the following function is a class function
    bool operator==(const DNA* & rhs){
        // since the uniqueness of an object is defined by sequence and location ID
        // the equality operator considers only those two criteria
        return ((getSequence() == rhs->getSequence()) && (getLocId() == rhs->getLocId()));
    }
    const DNA& operator=(const DNA& rhs){
        if (this != &rhs){
            m_sequence = rhs.m_sequence;
            m_location = rhs.m_location;
            m_used = rhs.m_used;
        }
        return *this;
    }
    private:
    string m_sequence;  // this is the object key
    int m_location;     // location ID that the DNA is found
    // the following variable is used for lazy delete scheme in hash table
    // if it is set to false, it means the bucket in the hash table is free for insert
    // if it is set to true, it means the bucket contains live data, and we cannot overwrite it
    bool m_used;
};
class DnaDb{
    public:
    friend class Grader;
    friend class Tester;
    DnaDb(int size, hash_fn hash, prob_t probing);
    ~DnaDb();
    // Returns Load factor of the new table
    float lambda() const;
    // Returns the ratio of deleted buckets in the new table
    float deletedRatio() const;
    // insert only happens in the new table
    bool insert(DNA dna);
    // remove can happen from either table
    bool remove(DNA dna);
    // find can happen in either table
    const DNA getDNA(string sequence, int location) const;
    // update the information
    bool updateLocId(DNA dna, int location);
    void changeProbPolicy(prob_t policy);
    void dump() const;
    private:
    hash_fn    m_hash;          // hash function
    prob_t     m_newPolicy;     // stores the change of policy request

    DNA**       m_currentTable; // hash table
    int        m_currentCap;    // hash table size (capacity)
    int        m_currentSize;   // current number of entries
                                // m_currentSize includes deleted entries 
    int        m_currNumDeleted;// number of deleted entries
    prob_t     m_currProbing;   // collision handling policy

    DNA**      m_oldTable;      // hash table
    int        m_oldCap;        // hash table size (capacity)
    int        m_oldSize;       // current number of entries
                                // m_oldSize includes deleted entries
    int        m_oldNumDeleted; // number of deleted entries
    prob_t     m_oldProbing;    // collision handling policy

    int        m_transferIndex; // this can be used as a temporary place holder
                                // during incremental transfer to scanning the table

    //private helper functions
    bool isPrime(int number);
    int findNextPrime(int current);

    /******************************************
    * Private function declarations go here! *
    ******************************************/

};
#endif