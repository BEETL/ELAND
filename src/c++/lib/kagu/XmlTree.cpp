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
 ** @file XmlTree.cpp
 **
 ** @brief A poor man's XML parser and writer while implementing a property map.
 **
 ** @author Michael Stromberg
 **/

#include "common/Exceptions.hh"
#include "common/StringUtilities.hh"
#include "kagu/XmlTree.h"

using namespace std;
namespace cc = casava::common;
using cc::StringUtilities;

namespace casava {
namespace kagu {

// our XML regular expressions
boost::regex XmlTree::mAttributeRegex("(\\S+?)=\"(.+?)\"");
boost::regex XmlTree::mElementRegex("^\\s*<(\\S+)\\s*([^/>]*)([/]*)>");
boost::regex XmlTree::mEndElementRegex("^\\s*</[^>]+>");
boost::regex XmlTree::mTextRegex("^(.+?)<");
boost::regex XmlTree::mXmlDeclarationRegex("^<\\?xml.+?>");

// constructor
XmlTree::XmlTree(void)
    : mpHead(NULL)
{
    mpHead = new XmlEntry();
    mpHead->Value = "<?xml version='1.0' standalone='yes'?>";
}

// destructor
XmlTree::~XmlTree(void) {
    DeleteTree(mpHead);
}

// adds data to the property tree
void XmlTree::AddStr(const string& key, const string& s) {
    vector<string> keys;
    StringUtilities::Split(key, '.', keys);

    XmlEntry* pEntry = mpHead;

    vector<XmlEntry*>::const_iterator xeIter;
    vector<string>::const_iterator sIter;

    for(sIter = keys.begin(); sIter != keys.end(); ++sIter) {

        // look for our key
        bool foundKey = false;
        for(xeIter = pEntry->Children.begin(); xeIter != pEntry->Children.end(); ++xeIter) {
            if((*xeIter)->Name == *sIter) {
                pEntry = *xeIter;
                foundKey = true;
                break;
            }
        }

        // create the key
        if(!foundKey) {
            XmlEntry* pChild = new XmlEntry();
            pEntry->Children.push_back(pChild);
            pChild->pParent = pEntry;
            pChild->Name = *sIter;
            pEntry = pChild;
        }
    }

    // assign the value
    pEntry->Value = s;
}

// adds data to the property tree
void XmlTree::AddInt(const string& key, const int32_t num) {
    ostringstream sb;
    sb << num;
    AddStr(key, sb.str());
}

void XmlTree::AddUInt(const string& key, const uint32_t num) {
    ostringstream sb;
    sb << num;
    AddStr(key, sb.str());
}

// adds data to the property tree
void XmlTree::AddDbl(const string& key, const double num) {
    ostringstream sb;
    sb << num;
    AddStr(key, sb.str());
}

// deletes a specified tree
void XmlTree::DeleteTree(XmlEntry* pEntry) {
    if(!pEntry) return;
    vector<XmlEntry*>::const_iterator xmpIter;
    for(xmpIter = pEntry->Children.begin(); xmpIter != pEntry->Children.end(); ++xmpIter) DeleteTree(*xmpIter);
    delete pEntry;
}

// returns all of the elements corresponding to the key
void XmlTree::GetElements(const string& key, Entries_t& entries) {

    // clear the supplied vector
    entries.clear();

    vector<string> keys;
    StringUtilities::Split(key, '.', keys);

    XmlEntry* pEntry = mpHead;

    // traverse down to the parent of the desired key
    vector<XmlEntry*>::const_iterator xeIter;
    vector<string>::const_iterator sIter;

    for(sIter = keys.begin(); sIter != (keys.end() - 1); ++sIter) {

        // look for our key
        bool foundKey = false;
        for(xeIter = pEntry->Children.begin(); xeIter != pEntry->Children.end(); ++xeIter) {
            if((*xeIter)->Name == *sIter) {
                pEntry = *xeIter;
                foundKey = true;
                break;
            }
        }

        // one of the keys leading up to the final key was missing
        if(!foundKey) return;
    }

    // create the entry vector
    for(xeIter = pEntry->Children.begin(); xeIter != pEntry->Children.end(); ++xeIter) {
        if((*xeIter)->Name == *sIter) {
            Entry e;
            e.Name       = (*xeIter)->Name;
            e.Value      = (*xeIter)->Value;
            e.Attributes = (*xeIter)->Attributes;
            entries.push_back(e);
        }
    }
}

// imports data from the XML file
void XmlTree::Import(const string& filename) {

    ifstream in(filename.c_str());

    if(in.fail()) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to open the XML file (%s) for reading") % filename).str()));
    }

    // initialization
    string line;
    queue<XmlNode> bufferQueue;
    XmlEntry* pCurrentEntry = mpHead;

    while(true) {

        if(bufferQueue.empty()) {

            // grab the next line
            getline(in, line);
            if(in.eof()) break;

            while(!line.empty()) {

                boost::smatch results;
                XmlNode node;
                uint32_t deleteLength = 0;

                if(boost::regex_search(line, results, mXmlDeclarationRegex)) {

                    // ===============
                    // XML declaration
                    // ===============

                    deleteLength  = (uint32_t)results[0].length();
                    node.Name     = "XmlDeclaration";
                    node.Value    = line.substr(0, deleteLength);
                    node.NodeType = XmlNodeType_XmlDeclaration;

                } else if(boost::regex_search(line, results, mEndElementRegex)) {

                    // ===========
                    // End Element
                    // ===========

                    deleteLength  = (uint32_t)results[0].length();
                    node.Name     = results[1].str();
                    node.NodeType = XmlNodeType_EndElement;

                } else if(boost::regex_search(line, results, mElementRegex)) {

                    // =======
                    // Element
                    // =======

                    deleteLength  = (uint32_t)results[0].length();
                    node.Name     = results[1].str();
                    node.NodeType = XmlNodeType_Element;

                    const bool hasEndElement = ((results.size() == 4) && (results[3] == "/"));

                    // extract attributes
                    const string& attribString = results[2].str();
                    boost::sregex_iterator attribIter(attribString.begin(), attribString.end(), mAttributeRegex);
                    boost::sregex_iterator attribEnd;

                    for(; attribIter != attribEnd; ++attribIter) {
                        KeyValue kv;
                        kv.Name  = (*attribIter)[1].str();
                        kv.Value = (*attribIter)[2].str();
                        node.Attributes.push_back(kv);
                    }

                    // handle an embedded "end element"
                    if(hasEndElement) {
                        bufferQueue.push(node);
                        node.Name     = results[1].str();
                        node.NodeType = XmlNodeType_EndElement;
                    }

                } else if(boost::regex_search(line, results, mTextRegex)) {

                    // ====
                    // Text
                    // ====

                    deleteLength  = (uint32_t)results[0].length() - 1;
                    node.Name     = "Text";
                    node.Value    = results[1].str();
                    node.NodeType = XmlNodeType_Text;

                } else {

                    // none of our regular expressions worked
                    BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("None of the XML regular expressions captured the current line: [%s]") % line).str()));
                }

                // add the node to the queue and remove the tag from the string
                bufferQueue.push(node);
                line.erase(0, deleteLength);
            }
        }

        // create the tree
        while(!bufferQueue.empty()) {

            XmlNode currentNode = bufferQueue.front();

            if(currentNode.NodeType == XmlNodeType_Element) {
                XmlEntry* pChild = new XmlEntry();
                pCurrentEntry->Children.push_back(pChild);
                pChild->pParent = pCurrentEntry;
                pCurrentEntry = pChild;
                pCurrentEntry->Name = currentNode.Name;
                if(!currentNode.Attributes.empty()) pCurrentEntry->Attributes = currentNode.Attributes;
            } else if(currentNode.NodeType == XmlNodeType_EndElement) {
                pCurrentEntry = pCurrentEntry->pParent;
            } else if(currentNode.NodeType == XmlNodeType_Text) {
                pCurrentEntry->Value = currentNode.Value;
            } else if(currentNode.NodeType == XmlNodeType_XmlDeclaration) {
                mpHead->Value = currentNode.Value;
            }

            bufferQueue.pop();
        }
    }

    in.close();
}

// dumps the data from the XML tree
void XmlTree::Print(void) {
    PrintTree(cout);
}

// prints from a tree node
void XmlTree::PrintBranch(ostream& out, const XmlEntry* pHead, uint32_t level) {

    string spacer(level * OUTPUT_INDENT_LEN, ' ');

    out << spacer << "<" << pHead->Name;
    if(!pHead->Attributes.empty()) {
        Attributes_t::const_iterator attribIter;
        for(attribIter = pHead->Attributes.begin(); attribIter != pHead->Attributes.end(); ++attribIter) {
            out << " " << attribIter->Name << "=\"" << attribIter->Value << "\"";
        }
    }
    out << ">";

    if(!pHead->Value.empty()) {
        out << pHead->Value << "</" << pHead->Name << ">" << endl;
    } else if(!pHead->Children.empty()) {
        out << endl;
    }

    if(!pHead->Children.empty()) {
        vector<XmlEntry*>::const_iterator childIter;
        for(childIter = pHead->Children.begin(); childIter != pHead->Children.end(); ++childIter) PrintBranch(out, *childIter, level + 1);
    }

    if(pHead->Value.empty()) {
        if(!pHead->Children.empty()) out << spacer;
        out << "</" << pHead->Name << ">" << endl;
    }
}

// prints the tree
void XmlTree::PrintTree(ostream& out) {
    if(!mpHead->Value.empty()) out << mpHead->Value << endl;
    vector<XmlEntry*>::const_iterator childIter;
    for(childIter = mpHead->Children.begin(); childIter != mpHead->Children.end(); ++childIter) {
        PrintBranch(out, *childIter, 0);
    }
}

// writes the content of the tree to a XML file
void XmlTree::Write(const string& filename) {
    ofstream out(filename.c_str(), ios::binary);

    if(out.fail()) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to open the XML file (%s) for writing") % filename).str()));
    }

    PrintTree(out);
    out.close();
}

}
}
