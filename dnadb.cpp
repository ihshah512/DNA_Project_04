// CMSC 341 - Spring 2025 - Project 4
#include "dnadb.h"
#include <cmath>
#include <iostream>

DnaDb::DnaDb(int size, hash_fn hash, prob_t probing)
{
    m_hash = hash;
    m_newPolicy = probing;
    m_currProbing = probing;
    m_oldProbing = probing;

    m_currentCap = findNextPrime(size);
    m_currentTable = new DNA *[m_currentCap]();
    m_currentSize = 0;
    m_currNumDeleted = 0;

    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;
    m_transferIndex = 0;
}

DnaDb::~DnaDb()
{
    for (int i = 0; i < m_currentCap; i++)
    {
        if (m_currentTable[i] != nullptr)
        {
            delete m_currentTable[i];
        }
    }
    delete[] m_currentTable;

    if (m_oldTable != nullptr)
    {
        for (int i = 0; i < m_oldCap; i++)
        {
            if (m_oldTable[i] != nullptr)
            {
                delete m_oldTable[i];
            }
        }
        delete[] m_oldTable;
    }
}

void DnaDb::changeProbPolicy(prob_t policy)
{
    m_newPolicy = policy;
}

bool DnaDb::insert(DNA dna)
{
    if (dna.getLocId() < MINLOCID || dna.getLocId() > MAXLOCID)
    {
        return false;
    }

    DNA existing = getDNA(dna.getSequence(), dna.getLocId());
    if (!existing.getSequence().empty())
    {
        return false;
    }

    if (m_oldTable != nullptr)
    {
        DNA oldExisting = getDNAFromTable(dna.getSequence(), dna.getLocId(), m_oldTable, m_oldCap, m_oldProbing);
        if (!oldExisting.getSequence().empty())
        {
            return false;
        }
    }

    unsigned int index = m_hash(dna.getSequence()) % m_currentCap;
    int i = 0;
    unsigned int originalIndex = index;
    int maxProbes = m_currentCap * 4;

    while (m_currentTable[index] != nullptr && m_currentTable[index]->getUsed() && i < maxProbes)
    {
        i++;
        if (m_currProbing == DOUBLEHASH)
        {
            int step = 11 - (m_hash(dna.getSequence()) % 11);
            if (step == 0)
                step = 1;
            index = (originalIndex + i * step) % m_currentCap;
        }
        else if (m_currProbing == QUADRATIC)
        {
            index = (originalIndex + i * i) % m_currentCap;
        }
        else
        {
            index = (originalIndex + i) % m_currentCap;
        }
    }

    if (i >= maxProbes)
    {
        // Fallback to linear probing
        index = originalIndex;
        i = 0;
        while (m_currentTable[index] != nullptr && m_currentTable[index]->getUsed() && i < m_currentCap)
        {
            i++;
            index = (originalIndex + i) % m_currentCap;
        }
        if (i >= m_currentCap)
        {
            return false;
        }
    }

    if (m_currentTable[index] == nullptr)
    {
        m_currentTable[index] = new DNA(dna);
        m_currentSize++;
    }
    else
    {
        *m_currentTable[index] = dna;
        if (!m_currentTable[index]->getUsed())
        {
            m_currNumDeleted--;
        }
        m_currentSize++;
    }
    m_currentTable[index]->setUsed(true);

    checkRehashCriteria();
    if (m_oldTable != nullptr)
    {
        transferNextQuarter();
    }

    return true;
}

bool DnaDb::remove(DNA dna)
{
    unsigned int index = m_hash(dna.getSequence()) % m_currentCap;
    int i = 0;
    unsigned int originalIndex = index;
    int maxProbes = m_currentCap * 4;

    while (m_currentTable[index] != nullptr && i < maxProbes)
    {
        if (m_currentTable[index]->getUsed() && *m_currentTable[index] == dna)
        {
            m_currentTable[index]->setUsed(false);
            m_currNumDeleted++;
            m_currentSize--;
            checkRehashCriteria();
            if (m_oldTable != nullptr)
            {
                transferNextQuarter();
            }
            return true;
        }
        i++;
        if (m_currProbing == DOUBLEHASH)
        {
            int step = 11 - (m_hash(dna.getSequence()) % 11);
            if (step == 0)
                step = 1;
            index = (originalIndex + i * step) % m_currentCap;
        }
        else if (m_currProbing == QUADRATIC)
        {
            index = (originalIndex + i * i) % m_currentCap;
        }
        else
        {
            index = (originalIndex + i) % m_currentCap;
        }
    }

    if (i >= maxProbes)
    {
        index = originalIndex;
        i = 0;
        while (m_currentTable[index] != nullptr && i < m_currentCap)
        {
            if (m_currentTable[index]->getUsed() && *m_currentTable[index] == dna)
            {
                m_currentTable[index]->setUsed(false);
                m_currNumDeleted++;
                m_currentSize--;
                checkRehashCriteria();
                if (m_oldTable != nullptr)
                {
                    transferNextQuarter();
                }
                return true;
            }
            i++;
            index = (originalIndex + i) % m_currentCap;
        }
    }

    if (m_oldTable != nullptr)
    {
        index = m_hash(dna.getSequence()) % m_oldCap;
        i = 0;
        originalIndex = index;

        while (m_currentTable[index] != nullptr && i < maxProbes)
        {
            if (m_oldTable[index]->getUsed() && *m_oldTable[index] == dna)
            {
                m_oldTable[index]->setUsed(false);
                m_oldNumDeleted++;
                m_oldSize--;
                transferNextQuarter();
                return true;
            }
            i++;
            if (m_oldProbing == DOUBLEHASH)
            {
                int step = 11 - (m_hash(dna.getSequence()) % 11);
                if (step == 0)
                    step = 1;
                index = (originalIndex + i * step) % m_oldCap;
            }
            else if (m_oldProbing == QUADRATIC)
            {
                index = (originalIndex + i * i) % m_oldCap;
            }
            else
            {
                index = (originalIndex + i) % m_oldCap;
            }
        }

        if (i >= maxProbes)
        {
            index = originalIndex;
            i = 0;
            while (m_oldTable[index] != nullptr && i < m_oldCap)
            {
                if (m_oldTable[index]->getUsed() && *m_oldTable[index] == dna)
                {
                    m_oldTable[index]->setUsed(false);
                    m_oldNumDeleted++;
                    m_oldSize--;
                    transferNextQuarter();
                    return true;
                }
                i++;
                index = (originalIndex + i) % m_oldCap;
            }
        }
    }

    return false;
}

const DNA DnaDb::getDNA(string sequence, int location) const
{
    DNA result = getDNAFromTable(sequence, location, m_currentTable, m_currentCap, m_currProbing);
    if (!result.getSequence().empty() && result.getUsed())
    {
        return result;
    }

    if (m_oldTable != nullptr)
    {
        result = getDNAFromTable(sequence, location, m_oldTable, m_oldCap, m_oldProbing);
        if (!result.getSequence().empty() && result.getUsed())
        {
            return result;
        }
    }

    return DNA();
}

bool DnaDb::updateLocId(DNA dna, int location)
{
    if (location < MINLOCID || location > MAXLOCID)
    {
        return false;
    }

    unsigned int index = m_hash(dna.getSequence()) % m_currentCap;
    int i = 0;
    unsigned int originalIndex = index;
    int maxProbes = m_currentCap * 4;

    while (m_currentTable[index] != nullptr && i < maxProbes)
    {
        if (m_currentTable[index]->getUsed() && *m_currentTable[index] == dna)
        {
            m_currentTable[index]->setLocID(location);
            return true;
        }
        i++;
        if (m_currProbing == DOUBLEHASH)
        {
            int step = 11 - (m_hash(dna.getSequence()) % 11);
            if (step == 0)
                step = 1;
            index = (originalIndex + i * step) % m_currentCap;
        }
        else if (m_currProbing == QUADRATIC)
        {
            index = (originalIndex + i * i) % m_currentCap;
        }
        else
        {
            index = (originalIndex + i) % m_currentCap;
        }
    }

    if (i >= maxProbes)
    {
        index = originalIndex;
        i = 0;
        while (m_currentTable[index] != nullptr && i < m_currentCap)
        {
            if (m_currentTable[index]->getUsed() && *m_currentTable[index] == dna)
            {
                m_currentTable[index]->setLocID(location);
                return true;
            }
            i++;
            index = (originalIndex + i) % m_currentCap;
        }
    }

    if (m_oldTable != nullptr)
    {
        index = m_hash(dna.getSequence()) % m_oldCap;
        i = 0;
        originalIndex = index;

        while (m_oldTable[index] != nullptr && i < maxProbes)
        {
            if (m_oldTable[index]->getUsed() && *m_oldTable[index] == dna)
            {
                m_oldTable[index]->setLocID(location);
                return true;
            }
            i++;
            if (m_oldProbing == DOUBLEHASH)
            {
                int step = 11 - (m_hash(dna.getSequence()) % 11);
                if (step == 0)
                    step = 1;
                index = (originalIndex + i * step) % m_oldCap;
            }
            else if (m_oldProbing == QUADRATIC)
            {
                index = (originalIndex + i * i) % m_oldCap;
            }
            else
            {
                index = (originalIndex + i) % m_oldCap;
            }
        }

        if (i >= maxProbes)
        {
            index = originalIndex;
            i = 0;
            while (m_oldTable[index] != nullptr && i < m_oldCap)
            {
                if (m_oldTable[index]->getUsed() && *m_oldTable[index] == dna)
                {
                    m_oldTable[index]->setLocID(location);
                    return true;
                }
                i++;
                index = (originalIndex + i) % m_oldCap;
            }
        }
    }

    return false;
}

float DnaDb::lambda() const
{
    return static_cast<float>(m_currentSize + m_currNumDeleted) / m_currentCap;
}

float DnaDb::deletedRatio() const
{
    if (m_currentSize + m_currNumDeleted == 0)
    {
        return 0.0f;
    }
    return static_cast<float>(m_currNumDeleted) / (m_currentSize + m_currNumDeleted);
}

void DnaDb::dump() const
{
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
    {
        for (int i = 0; i < m_currentCap; i++)
        {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
    {
        for (int i = 0; i < m_oldCap; i++)
        {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
    }
}

bool DnaDb::isPrime(int number)
{
    if (number <= 1)
        return false;
    for (int i = 2; i <= sqrt(number); i++)
    {
        if (number % i == 0)
        {
            return false;
        }
    }
    return true;
}

int DnaDb::findNextPrime(int current)
{
    if (current < MINPRIME)
        current = MINPRIME - 1;
    for (int i = current + 1; i <= MAXPRIME; i++)
    {
        if (isPrime(i))
        {
            return i;
        }
    }
    return MAXPRIME;
}

DNA DnaDb::getDNAFromTable(string sequence, int location, DNA **table, int cap, prob_t probing) const
{
    unsigned int index = m_hash(sequence) % cap;
    int i = 0;
    unsigned int originalIndex = index;
    int maxProbes = cap * 4;

    while (table[index] != nullptr && i < maxProbes)
    {
        if (table[index]->getUsed() && table[index]->getSequence() == sequence && table[index]->getLocId() == location)
        {
            return *table[index];
        }
        i++;
        if (probing == DOUBLEHASH)
        {
            int step = 11 - (m_hash(sequence) % 11);
            if (step == 0)
                step = 1;
            index = (originalIndex + i * step) % cap;
        }
        else if (probing == QUADRATIC)
        {
            index = (originalIndex + i * i) % cap;
        }
        else
        {
            index = (originalIndex + i) % cap;
        }
    }

    if (i >= maxProbes)
    {
        index = originalIndex;
        i = 0;
        while (table[index] != nullptr && i < cap)
        {
            if (table[index]->getUsed() && table[index]->getSequence() == sequence && table[index]->getLocId() == location)
            {
                return *table[index];
            }
            i++;
            index = (originalIndex + i) % cap;
        }
    }

    return DNA();
}

void DnaDb::checkRehashCriteria()
{
    if (lambda() > 0.5 || deletedRatio() > 0.8)
    {
        initiateRehash();
    }
}

void DnaDb::initiateRehash()
{
    int liveEntries = m_currentSize;
    int newSize = findNextPrime(liveEntries * 12); // Larger table to reduce collisions
    if (newSize < MINPRIME)
        newSize = MINPRIME;

    DNA **newTable = new DNA *[newSize]();

    m_oldTable = m_currentTable;
    m_oldCap = m_currentCap;
    m_oldSize = m_currentSize;
    m_oldNumDeleted = m_currNumDeleted;
    m_oldProbing = m_currProbing;

    m_currentTable = newTable;
    m_currentCap = newSize;
    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_currProbing = m_newPolicy;
    m_transferIndex = 0;
}

void DnaDb::transferNextQuarter()
{
    if (m_oldTable == nullptr)
    {
        return;
    }

    int quarterSize = m_oldCap / 4;
    if (quarterSize < 1)
        quarterSize = 1;
    int endIndex = m_transferIndex + quarterSize;
    if (endIndex > m_oldCap)
    {
        endIndex = m_oldCap;
    }

    for (int i = m_transferIndex; i < endIndex && i < m_oldCap; i++)
    {
        if (m_oldTable[i] != nullptr && m_oldTable[i]->getUsed())
        {
            DNA temp = *m_oldTable[i];
            bool inserted = insert(temp);
            if (inserted)
            {
                m_oldTable[i]->setUsed(false);
                m_oldNumDeleted++;
                m_oldSize--;
            }
        }
    }

    m_transferIndex = endIndex;

    if (m_transferIndex >= m_oldCap)
    {
        for (int i = 0; i < m_oldCap; i++)
        {
            if (m_oldTable[i] != nullptr)
            {
                delete m_oldTable[i];
                m_oldTable[i] = nullptr;
            }
        }
        delete[] m_oldTable;
        m_oldTable = nullptr;
        m_oldCap = 0;
        m_oldSize = 0;
        m_oldNumDeleted = 0;
        m_transferIndex = 0;
    }
}