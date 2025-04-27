// CMSC 341 - Spring 2025 - Project 4
#include "dnadb.h"
#include <math.h>
#include <algorithm>
#include <random>
#include <vector>
using namespace std;

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

enum RANDOM
{
    UNIFORMINT,
    UNIFORMREAL,
    NORMAL,
    SHUFFLE
};
class Random
{
public:
    Random() {}
    Random(int min, int max, RANDOM type = UNIFORMINT, int mean = 50, int stdev = 20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL)
        {
            // the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            // the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            // the mean and standard deviation can change by passing new values to constructor
            m_normdist = std::normal_distribution<>(mean, stdev);
        }
        else if (type == UNIFORMINT)
        {
            // the case of UNIFORMINT to generate integer numbers
            //  Using a fixed seed value generates always the same sequence
            //  of pseudorandom numbers, e.g. reproducing scientific experiments
            //  here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10); // 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min, max);
        }
        else if (type == UNIFORMREAL)
        {                                   // the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10); // 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min, (double)max);
        }
        else
        { // the case of SHUFFLE to generate every number only once
            m_generator = std::mt19937(m_device());
        }
    }
    void setSeed(int seedNum)
    {
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }
    void init(int min, int max)
    {
        m_min = min;
        m_max = max;
        m_type = UNIFORMINT;
        m_generator = std::mt19937(10); // 10 is the fixed seed value
        m_unidist = std::uniform_int_distribution<>(min, max);
    }
    void getShuffle(vector<int> &array)
    {
        // this function provides a list of all values between min and max
        // in a random order, this function guarantees the uniqueness
        // of every value in the list
        // the user program creates the vector param and passes here
        // here we populate the vector using m_min and m_max
        for (int i = m_min; i <= m_max; i++)
        {
            array.push_back(i);
        }
        shuffle(array.begin(), array.end(), m_generator);
    }

    void getShuffle(int array[])
    {
        // this function provides a list of all values between min and max
        // in a random order, this function guarantees the uniqueness
        // of every value in the list
        // the param array must be of the size (m_max-m_min+1)
        // the user program creates the array and pass it here
        vector<int> temp;
        for (int i = m_min; i <= m_max; i++)
        {
            temp.push_back(i);
        }
        std::shuffle(temp.begin(), temp.end(), m_generator);
        vector<int>::iterator it;
        int i = 0;
        for (it = temp.begin(); it != temp.end(); it++)
        {
            array[i] = *it;
            i++;
        }
    }

    int getRandNum()
    {
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if (m_type == NORMAL)
        {
            // returns a random number in a set with normal distribution
            // we limit random numbers by the min and max values
            result = m_min - 1;
            while (result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT)
        {
            // this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum()
    {
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result * 100.0) / 100.0;
        return result;
    }

    string getRandString(int size)
    {
        // the parameter size specifies the length of string we ask for
        // to use ASCII char the number range in constructor must be set to 97 - 122
        // and the Random type must be UNIFORMINT (it is default in constructor)
        string output = "";
        for (int i = 0; i < size; i++)
        {
            output = output + (char)getRandNum();
        }
        return output;
    }

    int getMin() { return m_min; }
    int getMax() { return m_max; }

private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;            // normal distribution
    std::uniform_int_distribution<> m_unidist;        // integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal; // real uniform distribution
};

class Tester
{
public:
    bool testInsertNonColliding()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);
        vector<DNA> inserted;
        int numInserts = 10;

        for (int i = 0; i < numInserts; i++)
        {
            DNA dna(sequencer(5, i), locGen.getRandNum(), true);
            inserted.push_back(dna);
            if (!db.insert(dna))
                return false;
            DNA found = db.getDNA(dna.getSequence(), dna.getLocId());
            if (found.getSequence() != dna.getSequence() ||
                found.getLocId() != dna.getLocId())
                return false;
        }

        float loadFactor = static_cast<float>(numInserts) / MINPRIME;
        if (fabs(db.lambda() - loadFactor) > 0.01f)
            return false;
        return true;
    }

    bool testGetDNANonExistent()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        DNA r = db.getDNA("AAAAA", MINLOCID);
        return r.getSequence().empty();
    }

    bool testGetDNANonColliding()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);
        vector<DNA> inserted;

        for (int i = 0; i < 5; i++)
        {
            DNA dna(sequencer(5, i), locGen.getRandNum(), true);
            inserted.push_back(dna);
            if (!db.insert(dna))
                return false;
        }
        for (auto &dna : inserted)
        {
            DNA f = db.getDNA(dna.getSequence(), dna.getLocId());
            if (f.getSequence() != dna.getSequence() ||
                f.getLocId() != dna.getLocId())
                return false;
        }
        return true;
    }

    bool testGetDNAColliding()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);
        vector<DNA> inserted;
        string seq = sequencer(5, 1);

        for (int i = 0; i < 3; i++)
        {
            DNA dna(seq, locGen.getRandNum(), true);
            inserted.push_back(dna);
            if (!db.insert(dna))
                return false;
        }
        for (auto &dna : inserted)
        {
            DNA f = db.getDNA(dna.getSequence(), dna.getLocId());
            if (f.getSequence() != dna.getSequence() ||
                f.getLocId() != dna.getLocId())
                return false;
        }
        return true;
    }

    bool testRemoveNonColliding()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);
        vector<DNA> inserted;

        for (int i = 0; i < 5; i++)
        {
            DNA dna(sequencer(5, i), locGen.getRandNum(), true);
            inserted.push_back(dna);
            if (!db.insert(dna))
                return false;
        }
        for (auto &dna : inserted)
        {
            if (!db.remove(dna))
                return false;
            DNA f = db.getDNA(dna.getSequence(), dna.getLocId());
            if (!f.getSequence().empty())
                return false;
        }
        return true;
    }

    bool testRemoveColliding()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);
        vector<DNA> inserted;
        string seq = sequencer(5, 1);

        for (int i = 0; i < 3; i++)
        {
            DNA dna(seq, locGen.getRandNum(), true);
            inserted.push_back(dna);
            if (!db.insert(dna))
                return false;
        }
        for (auto &dna : inserted)
        {
            if (!db.remove(dna))
                return false;
            DNA f = db.getDNA(dna.getSequence(), dna.getLocId());
            if (!f.getSequence().empty())
                return false;
        }
        return true;
    }

    bool testRehashLoadFactor()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);

        // trigger at just over 50%
        int threshold = static_cast<int>(MINPRIME * 0.5) + 1;
        for (int i = 0; i < threshold; i++)
        {
            DNA dna(sequencer(5, i), locGen.getRandNum(), true);
            if (!db.insert(dna))
                return false;
        }
        if (db.m_oldTable == nullptr)
            return false;

        float loadFactor = static_cast<float>(threshold) / MINPRIME;
        if (fabs(db.lambda() - loadFactor) > 0.01f)
            return false;
        return true;
    }

    bool testRehashDeleteRatio()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);
        vector<DNA> inserted;

        int numInserts = 20;
        for (int i = 0; i < numInserts; i++)
        {
            DNA dna(sequencer(5, i), locGen.getRandNum(), true);
            inserted.push_back(dna);
            if (!db.insert(dna))
                return false;
        }

        int numToDelete = static_cast<int>(numInserts * 0.8);
        for (int i = 0; i < numToDelete; i++)
        {
            if (!db.remove(inserted[i]))
                return false;
        }

        if (db.m_oldTable == nullptr)
            return false;
        return true;
    }
    bool testRehashCompletionDeleteRatio()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);

        int numInserts = 20;
        vector<DNA> inserted;
        for (int i = 0; i < numInserts; ++i)
        {
            DNA dna(sequencer(5, i), locGen.getRandNum(), true);
            inserted.push_back(dna);
            if (!db.insert(dna))
                return false;
        }
        int numToDelete = static_cast<int>(numInserts * 0.8);
        for (int i = 0; i < numToDelete; ++i)
        {
            if (!db.remove(inserted[i]))
                return false;
        }
        if (!db.m_oldTable)
            return false;

        // Keep inserting until old table is fully transferred
        int extra = 0;
        while (db.m_oldTable)
        {
            DNA dummy(sequencer(5, numInserts + extra), locGen.getRandNum(), true);
            if (!db.insert(dummy))
                return false;
            ++extra;
        }
        if (db.m_oldTable)
            return false;

        // Verify remaining originals still present
        for (int i = numToDelete; i < numInserts; ++i)
        {
            DNA f = db.getDNA(inserted[i].getSequence(), inserted[i].getLocId());
            if (f.getSequence() != inserted[i].getSequence() ||
                f.getLocId() != inserted[i].getLocId())
                return false;
        }
        return true;
    }
    bool testRehashCompletionLoadFactor()
    {
        DnaDb db(MINPRIME, hashCode, DOUBLEHASH);
        Random locGen(MINLOCID, MAXLOCID);
        locGen.setSeed(42);

        int threshold = static_cast<int>(MINPRIME * 0.5) + 1;
        vector<DNA> inserted;
        for (int i = 0; i < threshold; ++i)
        {
            DNA dna(sequencer(5, i), locGen.getRandNum(), true);
            inserted.push_back(dna);
            if (!db.insert(dna))
                return false;
        }
        // Now rehash has started
        if (!db.m_oldTable)
            return false;

        // Keep inserting until old table is fully transferred
        int extra = 0;
        while (db.m_oldTable)
        {
            DNA dummy(sequencer(5, threshold + extra), locGen.getRandNum(), true);
            if (!db.insert(dummy))
                return false;
            ++extra;
        }
        // Verify nothing left in old table
        if (db.m_oldTable)
            return false;

        // All original entries still found
        for (auto &dna : inserted)
        {
            DNA f = db.getDNA(dna.getSequence(), dna.getLocId());
            if (f.getSequence() != dna.getSequence() ||
                f.getLocId() != dna.getLocId())
                return false;
        }
        return true;
    }
};

int main()
{
    vector<DNA> dataList;
    Random RndLocation(MINLOCID, MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH);
    bool result = true;

    cout << "Inserting 49 data nodes!" << endl;
    for (int i = 0; i < 49; i++)
    {
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        // saving data for later use
        dataList.push_back(dataObj);
        // inserting data in to the DB object
        if (!dnadb.insert(dataObj))
            cout << "Did not insert " << &dataObj << endl;
    }
    // try to delete some data node
    cout << "Removing data node " << &dataList[5] << endl;
    dnadb.remove(dataList[5]);
    cout << "Removing data node " << &dataList[15] << endl;
    dnadb.remove(dataList[15]);
    dnadb.dump();

    // checking whether all data points are there
    cout << endl
         << "Checking whether all data exist in the DB:" << endl;
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++)
    {
        DNA anObj = dnadb.getDNA((*it).getSequence(), (*it).getLocId());
        bool foundIt = (*it == anObj);
        result = result && foundIt;
        if (!foundIt)
        {
            cout << "Data point " << (*it).getSequence() << "(" << (*it).getLocId() << ")" << " is missing!" << endl;
        }
    }
    if (result)
        cout << "\tAll data points exist in the DnaDb object!\n";

    Tester tester;
    cout << "Testing DnaDb implementation:\n";
    cout << "testInsertNonColliding: " << (tester.testInsertNonColliding() ? "PASS\n" : "FAIL\n");
    cout << "testGetDNANonExistent: " << (tester.testGetDNANonExistent() ? "PASS\n" : "FAIL\n");
    cout << "testGetDNANonColliding: " << (tester.testGetDNANonColliding() ? "PASS\n" : "FAIL\n");
    cout << "testGetDNAColliding: " << (tester.testGetDNAColliding() ? "PASS\n" : "FAIL\n");
    cout << "testRemoveNonColliding: " << (tester.testRemoveNonColliding() ? "PASS\n" : "FAIL\n");
    cout << "testRemoveColliding: " << (tester.testRemoveColliding() ? "PASS\n" : "FAIL\n");
    cout << "testRehashLoadFactor: " << (tester.testRehashLoadFactor() ? "PASS\n" : "FAIL\n");
    cout << "testRehashCompletionLoadFactor: " << (tester.testRehashCompletionLoadFactor() ? "PASS\n" : "FAIL\n");
    cout << "testRehashDeleteRatio: " << (tester.testRehashDeleteRatio() ? "PASS\n" : "FAIL\n");
    cout << "testRehashCompletionDeleteRatio: " << (tester.testRehashCompletionDeleteRatio() ? "PASS\n" : "FAIL\n");

    return 0;
}

unsigned int hashCode(const string str)
{
    unsigned int val = 0;
    const unsigned int thirtyThree = 33; // magic number from textbook
    for (int i = 0; i < str.length(); i++)
        val = val * thirtyThree + str[i];
    return val;
}
string sequencer(int size, int seedNum)
{
    // this function returns a random DNA sequence
    //  size param specifies the size of string
    string sequence = "";
    Random rndObject(0, 3);
    rndObject.setSeed(seedNum);
    for (int i = 0; i < size; i++)
    {
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}