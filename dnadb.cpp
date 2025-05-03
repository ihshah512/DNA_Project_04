#include "dnadb.h"
#include <iostream> // Only used for dump(), remove if not allowed

// parametierized constructor
DnaDb::DnaDb(int size, hash_fn hash, prob_t probing = DEFPOLCY)
{
    m_hash = hash;           // type of this function is defined in header file
    m_currProbing = probing; // This would default which is quaradtic
    m_newPolicy = probing;   // stores changes to policy request
    m_oldTable = nullptr;
    m_oldCap = 0;  // old hash table capacity
    m_oldSize = 0; // current numbner of entires
    m_oldNumDeleted = 0;
    m_transferIndex = 0;
    m_oldProbing = probing; // collision handling policy

    // making sure the prime bounds are within prescribed limits
    if (size < MINPRIME)
        size = MINPRIME;
    if (size > MAXPRIME)
        size = MAXPRIME;
    if (!isPrime(size))             // in case user input is not a prime number
        size = findNextPrime(size); // usage of helper funcitn to obtain prime number

    m_currentCap = size;
    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_currentTable = new DNA *[m_currentCap]; // memory initialization for the table
    for (int i = 0; i < m_currentCap; i++)
        m_currentTable[i] = nullptr;
}

// Destructor
DnaDb::~DnaDb()
{
    // Clean up current table if it exists and is different from old table
    if (m_currentTable && m_currentTable != m_oldTable) // if curent table exists and its not equal to old table
    {
        for (int i = 0; i < m_currentCap; i++)
        {
            delete m_currentTable[i];
        }
        delete[] m_currentTable;
    }

    // Clean up old table if it exists and is different from current table
    if (m_oldTable)
    {
        for (int i = 0; i < m_oldCap; i++)
        {
            delete m_oldTable[i];
        }
        delete[] m_oldTable;
    }
}

bool DnaDb::insert(DNA dna)
{
    // either location ID (not in range minlocid and maxloc id) or DNA sequence is empty
    if (dna.getLocId() < MINLOCID || dna.getLocId() > MAXLOCID || dna.getSequence().empty())
        return false;
    // no duplicates
    if (getDNA(dna.getSequence(), dna.getLocId()).getSequence() != "")
        return false;

    // if mid-rehash, move next quarter first
    if (m_oldTable)
        transferNextQuarter();

    unsigned idx = m_hash(dna.getSequence()) % m_currentCap;
    unsigned i = 0, pos = idx;
    while (true)
    {
        if (!m_currentTable[pos] || !m_currentTable[pos]->m_used)
        {
            if (m_currentTable[pos])
            {
                delete m_currentTable[pos];
                m_currentTable[pos] = nullptr;
            }
            m_currentTable[pos] = new DNA(dna);
            m_currentTable[pos]->m_used = true; // marking the position used
            ++m_currentSize;
            break;
        }
        ++i;
        if (m_currProbing == LINEAR)
            pos = (idx + i) % m_currentCap;
        else if (m_currProbing == QUADRATIC)
            pos = (idx + i * i) % m_currentCap;
        else // DOUBLEHASH
            pos = (idx + i * (11 - (m_hash(dna.getSequence()) % 11))) % m_currentCap;
    }

    checkRehashCriteria();
    return true;
}

// Public remove
bool DnaDb::remove(DNA dna)
{
    if (m_oldTable)
        transferNextQuarter();

    // search current table
    {
        unsigned idx = m_hash(dna.getSequence()) % m_currentCap;
        unsigned i = 0, pos = idx;
        while (m_currentTable[pos])
        {
            if (m_currentTable[pos]->m_used && *m_currentTable[pos] == dna)
            {
                m_currentTable[pos]->m_used = false;
                ++m_currNumDeleted;
                checkRehashCriteria();
                return true;
            }
            ++i;
            if (m_currProbing == LINEAR)
                pos = (idx + i) % m_currentCap;
            else if (m_currProbing == QUADRATIC)
                pos = (idx + i * i) % m_currentCap;
            else
                pos = (idx + i * (11 - (m_hash(dna.getSequence()) % 11))) % m_currentCap;
        }
    }
    // search old table
    if (m_oldTable)
    {
        unsigned idx = m_hash(dna.getSequence()) % m_oldCap;
        unsigned i = 0, pos = idx;
        while (m_oldTable[pos])
        {
            if (m_oldTable[pos]->m_used && *m_oldTable[pos] == dna)
            {
                m_oldTable[pos]->m_used = false;
                ++m_oldNumDeleted;
                checkRehashCriteria();
                return true;
            }
            ++i;
            if (m_oldProbing == LINEAR)
                pos = (idx + i) % m_oldCap;
            else if (m_oldProbing == QUADRATIC)
                pos = (idx + i * i) % m_oldCap;
            else
                pos = (idx + i * (11 - (m_hash(dna.getSequence()) % 11))) % m_oldCap;
        }
    }
    return false;
}

// dna lookup with given location and sequence
const DNA DnaDb::getDNA(string sequence, int location) const
{
    DNA found = getDNAFromTable(sequence, location,
                                m_currentTable, m_currentCap, m_currProbing);
    if (found.getSequence().empty() && m_oldTable)
        found = getDNAFromTable(sequence, location,
                                m_oldTable, m_oldCap, m_oldProbing);
    return found;
}

//  update location id based on seq and id
bool DnaDb::updateLocId(DNA dna, int location)
{
    DNA oldDNA = getDNA(dna.getSequence(), dna.getLocId());
    if (oldDNA.getSequence().empty())
        return false;
    remove(dna);
    DNA updated = dna;
    updated.setLocID(location);
    insert(updated);
    return true;
}

// Load factor — during rehash show old-table’s load, otherwise current
float DnaDb::lambda() const
{
    if (m_oldTable)
        return static_cast<float>(m_oldSize) / m_oldCap;
    return static_cast<float>(m_currentSize) / m_currentCap;
}

// Deleted-ratio — same: old table if rehashing, else current
float DnaDb::deletedRatio() const
{
    if (m_oldTable)
    {
        if (m_oldSize == 0)
            return 0.0f;
        return static_cast<float>(m_oldNumDeleted) / m_oldSize;
    }
    if (m_currentSize == 0)
        return 0.0f;
    return static_cast<float>(m_currNumDeleted) / m_currentSize;
}

// change policy
void DnaDb::changeProbPolicy(prob_t policy)
{
    m_newPolicy = policy;
}

// debug dump
void DnaDb::dump() const
{
    cout << "Dump for the current table: " << endl;
    for (int i = 0; i < m_currentCap; ++i)
        cout << "[" << i << "] : " << m_currentTable[i] << "\n";
    if (m_oldTable)
    {
        cout << "Dump for the old table: " << endl;
        for (int i = 0; i < m_oldCap; ++i)
            cout << "[" << i << "] : " << m_oldTable[i] << "\n";
    }
}
/*********************************************************** *|
 *             HELPER FUNCTIONS ARE BELOW                     |
 *                                                            /
 *                                                            /
 * ***********************************************************/
// Checking a number is prime or not
bool DnaDb::isPrime(int number)
{
    if (number <= 1)
        return false;
    if (number == 2)
        return true;
    if (number % 2 == 0)
        return false;
    for (int i = 3; i * i <= number; i += 2)
        if (number % i == 0)
            return false;
    return true;
}

// Find smallest prime > current, bounded by MINPRIME/MAXPRIME
int DnaDb::findNextPrime(int current)
{
    if (current <= MINPRIME)
        return MINPRIME;
    if (current >= MAXPRIME)
        return MAXPRIME;
    int next = current;
    while (!isPrime(next))
        ++next;
    return next;
}

// table-scan helper this function will be used as healper function in getDNA function
DNA DnaDb::getDNAFromTable(string sequence, int location,
                           DNA **table, int cap, prob_t probing) const
{
    unsigned idx = m_hash(sequence) % cap;
    unsigned i = 0, pos = idx;
    while (table[pos])
    {
        if (table[pos]->m_used &&
            table[pos]->m_sequence == sequence &&
            table[pos]->m_location == location)
            return *table[pos];
        ++i;
        if (probing == LINEAR)
            pos = (idx + i) % cap;
        else if (probing == QUADRATIC)
            pos = (idx + i * i) % cap;
        else
            pos = (idx + i * (11 - (m_hash(sequence) % 11))) % cap;
    }
    return DNA();
}

// rehash criteria
void DnaDb::checkRehashCriteria()
{
    // only start a rehash if one isn’t already in progress
    if (!m_oldTable && (lambda() > 0.5f || deletedRatio() >= 0.8f))
    {
        initiateRehash();
    }
}

// start a new incremental rehash
void DnaDb::initiateRehash()
{
    m_oldTable = m_currentTable;
    m_oldCap = m_currentCap;
    m_oldSize = m_currentSize;
    m_oldNumDeleted = m_currNumDeleted;
    m_oldProbing = m_currProbing;

    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_transferIndex = 0;
    m_currentTable = nullptr;

    int newCap = findNextPrime((m_oldSize - m_oldNumDeleted) * 4);
    m_currentCap = newCap;
    m_currProbing = m_newPolicy;
    m_currentTable = new DNA *[m_currentCap];
    for (int i = 0; i < m_currentCap; i++)
        m_currentTable[i] = nullptr;
}

// helper used during transfer
void DnaDb::rehashInsert(const DNA &dna)
{
    // cout << "rehashInsert: sequence=" << dna.getSequence() << ", loc=" << dna.getLocId() << endl;
    unsigned idx = m_hash(dna.getSequence()) % m_currentCap;
    unsigned i = 0, pos = idx;
    unsigned max_attempts = m_currentCap; // Prevent infinite loop
    while (i < max_attempts)
    {
        if (!m_currentTable[pos] || !m_currentTable[pos]->m_used)
        {
            if (m_currentTable[pos])
            {
                // cout << "Deleting existing at pos=" << pos << endl;
                delete m_currentTable[pos];
                m_currentTable[pos] = nullptr;
            }
            m_currentTable[pos] = new DNA(dna);
            m_currentTable[pos]->m_used = true;
            ++m_currentSize;
            // cout << "Inserted at pos=" << pos << endl;
            return;
        }
        ++i;
        if (m_currProbing == LINEAR)
            pos = (idx + i) % m_currentCap;
        else if (m_currProbing == QUADRATIC)
            pos = (idx + i * i) % m_currentCap;
        else
            pos = (idx + i * (11 - (m_hash(dna.getSequence()) % 11))) % m_currentCap;
    }
    // Table full; return without inserting (should not happen due to rehash sizing)
    // cout << "Warning: Cannot insert, table full" << endl;
    return;
}

// move 25% of old‐table entries
void DnaDb::transferNextQuarter()
{
    if (!m_oldTable)
        return;
    int quarter = m_oldCap / 4; // 25 of old table
    int start = m_transferIndex;
    int end = start + quarter < m_oldCap ? start + quarter : m_oldCap; // Replace std::min
    // cout << "Transferring from " << start << " to " << end << endl;
    for (int i = start; i < end; ++i)
    {
        if (m_oldTable[i] && m_oldTable[i]->m_used)
        {
            // cout << "Transferring DNA at index " << i << ": seq=" << m_oldTable[i]->getSequence() << endl;
            rehashInsert(*m_oldTable[i]);
            m_oldTable[i]->m_used = false;
            // Mark as transferred but don’t delete yet
        }
    }
    m_transferIndex = end;
    if (m_transferIndex >= m_oldCap)
    {
        // cout << "Cleaning up old table" << endl;
        for (int i = 0; i < m_oldCap; i++)
        {
            if (m_oldTable[i])
            {
                delete m_oldTable[i];
                m_oldTable[i] = nullptr;
            }
        }
        delete[] m_oldTable;
        m_oldTable = nullptr;
        m_oldCap = m_oldSize = m_oldNumDeleted = m_transferIndex = 0;
    }
}

// Copy Constructor
DnaDb::DnaDb(const DnaDb &other) : m_hash(other.m_hash),
                                   m_currProbing(other.m_currProbing),
                                   m_newPolicy(other.m_newPolicy),
                                   m_currentCap(other.m_currentCap),
                                   m_currentSize(other.m_currentSize),
                                   m_currNumDeleted(other.m_currNumDeleted),
                                   m_oldCap(other.m_oldCap),
                                   m_oldSize(other.m_oldSize),
                                   m_oldNumDeleted(other.m_oldNumDeleted),
                                   m_transferIndex(other.m_transferIndex),
                                   m_oldProbing(other.m_oldProbing),
                                   m_currentTable(nullptr),
                                   m_oldTable(nullptr)
{
    // Copy current table
    if (other.m_currentTable)
    {
        m_currentTable = new DNA *[m_currentCap];
        for (int i = 0; i < m_currentCap; i++)
        {
            if (other.m_currentTable[i] && other.m_currentTable[i]->m_used)
            {
                m_currentTable[i] = new DNA(*other.m_currentTable[i]);
            }
            else
            {
                m_currentTable[i] = nullptr;
            }
        }
    }

    // Copy old table if it exists
    if (other.m_oldTable)
    {
        m_oldTable = new DNA *[m_oldCap];
        for (int i = 0; i < m_oldCap; i++)
        {
            if (other.m_oldTable[i] && other.m_oldTable[i]->m_used)
            {
                m_oldTable[i] = new DNA(*other.m_oldTable[i]);
            }
            else
            {
                m_oldTable[i] = nullptr;
            }
        }
    }
}

// Assignment Operator (without swap)
DnaDb &DnaDb::operator=(const DnaDb &other)
{
    if (this != &other)
    { // Check for self-assignment
        // Clean up existing resources
        if (m_currentTable)
        {
            for (int i = 0; i < m_currentCap; i++)
            {
                delete m_currentTable[i];
            }
            delete[] m_currentTable;
        }
        if (m_oldTable)
        {
            for (int i = 0; i < m_oldCap; i++)
            {
                delete m_oldTable[i];
            }
            delete[] m_oldTable;
        }

        // Copy primitive members
        m_hash = other.m_hash;
        m_currProbing = other.m_currProbing;
        m_newPolicy = other.m_newPolicy;
        m_currentCap = other.m_currentCap;
        m_currentSize = other.m_currentSize;
        m_currNumDeleted = other.m_currNumDeleted;
        m_oldCap = other.m_oldCap;
        m_oldSize = other.m_oldSize;
        m_oldNumDeleted = other.m_oldNumDeleted;
        m_transferIndex = other.m_transferIndex;
        m_oldProbing = other.m_oldProbing;

        // Copy current table
        m_currentTable = new DNA *[m_currentCap];
        for (int i = 0; i < m_currentCap; i++)
        {
            m_currentTable[i] = nullptr;
            if (other.m_currentTable[i] && other.m_currentTable[i]->m_used)
            {
                m_currentTable[i] = new DNA(*other.m_currentTable[i]);
            }
        }

        // Copy old table if it exists
        m_oldTable = nullptr;
        if (other.m_oldTable)
        {
            m_oldTable = new DNA *[m_oldCap];
            for (int i = 0; i < m_oldCap; i++)
            {
                m_oldTable[i] = nullptr;
                if (other.m_oldTable[i] && other.m_oldTable[i]->m_used)
                {
                    m_oldTable[i] = new DNA(*other.m_oldTable[i]);
                }
            }
        }
    }
    return *this;
}
