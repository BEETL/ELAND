/**
 ** Copyright (c) 2007-2010 Illumina, Inc.
 **
 ** This software is covered by the "Illumina Genome Analyzer Software
 ** License Agreement" and the "Illumina Source Code License Agreement",
 ** and certain third party copyright/licenses, and any user of this
 ** source file is bound by the terms therein (see accompanying files
 ** Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
 ** Illumina_Source_Code_License_Agreement.pdf and third party
 ** copyright/license notices).
 **
 ** This file is part of the Consensus Assessment of Sequence And VAriation
 ** (CASAVA) software package.
 **
 ** @file XmlTree.h
 **
 ** @brief A poor man's XML parser and writer while implementing a property map.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <boost/algorithm/string/classification.hpp>
#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <stdint.h>
#include <string>

#define OUTPUT_INDENT_LEN 2

namespace casava {
namespace kagu {

enum XmlNodeType {
    XmlNodeType_Element,
    XmlNodeType_EndElement,
    XmlNodeType_Text,
    XmlNodeType_XmlDeclaration
};

struct KeyValue {
    std::string Name;
    std::string Value;
};

typedef std::vector<KeyValue> Attributes_t;

struct Entry : KeyValue {
    Attributes_t Attributes;
};

typedef std::vector<Entry> Entries_t;

struct XmlNode : Entry {
    XmlNodeType NodeType;
};

struct XmlEntry : Entry {
    XmlEntry* pParent;
    std::vector<XmlEntry*> Children;

    // constructor
    XmlEntry(void)
        : pParent(NULL)
    {}
};

class XmlTree {
public:
    // constructor
    XmlTree(void);
    // destructor
    ~XmlTree(void);
    // adds data to the property tree
    void AddStr(const std::string& key, const std::string& s);
    void AddInt(const std::string& key, const int32_t num);
    void AddUInt(const std::string& key, const uint32_t num);
    void AddDbl(const std::string& key, const double num);
    // returns all of the elements corresponding to the key
    void GetElements(const std::string& key, Entries_t& entries);
    // imports data from the XML file
    void Import(const std::string& filename);
    // dumps the contents of the tree to stdout
    void Print(void);
    // writes the contents of the tree to a XML file
    void Write(const std::string& filename);

private:
    // deletes a specified tree
    void DeleteTree(XmlEntry* pHead);
    // prints the tree starting at the specified branch
    void PrintBranch(std::ostream& out, const XmlEntry* pHead, uint32_t level);
    // prints the tree
    void PrintTree(std::ostream& out);
    // our head entry
    XmlEntry* mpHead;
    // our XML regular expressions
    static boost::regex mAttributeRegex;
    static boost::regex mElementRegex;
    static boost::regex mEndElementRegex;
    static boost::regex mTextRegex;
    static boost::regex mXmlDeclarationRegex;
};

}
}
