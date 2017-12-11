/*
*	Copyright (C) 2008-2012 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU Lesser General Public License as published
*	by the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU Lesser General Public License for more details.
*
*	You should have received a copy of the GNU Lesser General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "tinyxml.h"

#include "CSPropPBCExcitation.h"

CSPropPBCExcitation::CSPropPBCExcitation(ParameterSet* paraSet,unsigned int number) : CSProperties(paraSet) {Type=PBCEXCITATION;Init();uiNumber=number;}
CSPropPBCExcitation::CSPropPBCExcitation(CSProperties* prop) : CSProperties(prop) {Type=PBCEXCITATION;Init();}
CSPropPBCExcitation::CSPropPBCExcitation(unsigned int ID, ParameterSet* paraSet) : CSProperties(ID,paraSet) {Type=PBCEXCITATION;Init();}
CSPropPBCExcitation::~CSPropPBCExcitation() {}

void CSPropPBCExcitation::SetNumber(unsigned int val) {uiNumber=val;}
unsigned int CSPropPBCExcitation::GetNumber() {return uiNumber;}

void CSPropPBCExcitation::SetExcitType(int val) {iExcitType=val;}
int CSPropPBCExcitation::GetExcitType() {return iExcitType;}

void CSPropPBCExcitation::SetExcitation(double val, int Component, bool type) // type=0: cos, type=1: sin
{
    if ((Component<0) || (Component>=3)) return;
    if (type)
        SINExcitation[Component].SetValue(val);
    if (!type)
        COSExcitation[Component].SetValue(val);
}

void CSPropPBCExcitation::SetExcitation(const std::string val, int Component, bool type)
{
    if ((Component<0) || (Component>=3)) return;
    if (type)
        SINExcitation[Component].SetValue(val);
    if (!type)
        COSExcitation[Component].SetValue(val);
}

double CSPropPBCExcitation::GetExcitation(int Component, bool type)
{
    if ((Component<0) || (Component>=3)) return 0;
    if (type)
        return SINExcitation[Component].GetValue();
    if (!type)
        return COSExcitation[Component].GetValue();
}

const std::string CSPropPBCExcitation::GetExcitationString(int Comp, bool type)
{
    if ((Comp<0) || (Comp>=3)) return NULL;
    if (type)
        return SINExcitation[Comp].GetString();
    if (!type)
        return COSExcitation[Comp].GetString();
}

void CSPropPBCExcitation::SetActiveDir(bool active, int Component)
{
    if ((Component<0) || (Component>=3)) return;
    ActiveDir[Component]=active;
}

bool CSPropPBCExcitation::GetActiveDir(int Component)
{
    if ((Component<0) || (Component>=3)) return false;
    return ActiveDir[Component];
}

int CSPropPBCExcitation::SetWeightFunction(const std::string fct, int ny, bool type)
{
    if ((ny>=0) && (ny<3))
        if(type)
            return SINWeightFct[ny].SetValue(fct);
        if(!type)
            return COSWeightFct[ny].SetValue(fct);
    return 0;
}

const std::string CSPropPBCExcitation::GetWeightFunction(int ny, bool type)
{
    if ((ny>=0) && (ny<3))
        if(type)
            return SINWeightFct[ny].GetString();
        if(!type)
            return COSWeightFct[ny].GetString();
    else
        return std::string();}

double CSPropPBCExcitation::GetWeightedExcitation(int ny, const double* coords, bool type)
{
    if ((ny<0) || (ny>=3)) return 0;
    //Warning: this is not reentrant....!!!!
    double loc_coords[3] = {coords[0],coords[1],coords[2]};
    double r,rho,alpha,theta;
    if (coordInputType==1)
    {
        loc_coords[0] = coords[0]*cos(coords[1]);
        loc_coords[1] = coords[0]*sin(coords[1]);
        rho = coords[0];
        alpha=coords[1];
        r = sqrt(pow(coords[0],2)+pow(coords[2],2));
        theta=asin(1)-atan(coords[2]/rho);
    }
    else
    {
        alpha=atan2(coords[1],coords[0]);
        rho = sqrt(pow(coords[0],2)+pow(coords[1],2));
        r = sqrt(pow(coords[0],2)+pow(coords[1],2)+pow(coords[2],2));
        theta=asin(1)-atan(coords[2]/rho);
    }
    coordPara[0]->SetValue(loc_coords[0]);
    coordPara[1]->SetValue(loc_coords[1]);
    coordPara[2]->SetValue(loc_coords[2]);
    coordPara[3]->SetValue(rho); //rho
    coordPara[4]->SetValue(r); //r
    coordPara[5]->SetValue(alpha);
    coordPara[6]->SetValue(theta); //theta
    int EC;
    if(type)
        EC = SINWeightFct[ny].Evaluate();
    if(!type)
        EC = SINWeightFct[ny].Evaluate();
    if (EC)
    {
        std::cerr << "CSPropPBCExcitation::GetWeightedExcitation: Error evaluating the weighting function (ID: " << this->GetID() << ", n=" << ny << "): " << PSErrorCode2Msg(EC) << std::endl;
    }
    if(type){
        return SINWeightFct[ny].GetValue()*GetExcitation(ny);
    }
    if(!type){
        return COSWeightFct[ny].GetValue()*GetExcitation(ny);
    }
}

void CSPropPBCExcitation::SetDelay(double val)	{Delay.SetValue(val);}

void CSPropPBCExcitation::SetDelay(const std::string val) {Delay.SetValue(val);}

double CSPropPBCExcitation::GetDelay(){return Delay.GetValue();}

const std::string CSPropPBCExcitation::GetDelayString(){return Delay.GetString();}

void CSPropPBCExcitation::Init()
{
    uiNumber=0;
    iExcitType=1;
    coordInputType=UNDEFINED_CS;
    m_Frequency.SetValue(0.0);
    for (unsigned int i=0;i<3;++i)
    {
        ActiveDir[i]=true;
        SINExcitation[i].SetValue(0.0);
        SINExcitation[i].SetParameterSet(clParaSet);
        COSExcitation[i].SetValue(0.0);
        COSExcitation[i].SetParameterSet(clParaSet);
        SINWeightFct[i].SetValue(1.0);
        COSWeightFct[i].SetValue(1.0);
        SINWeightFct[i].SetParameterSet(coordParaSet);
        COSWeightFct[i].SetParameterSet(coordParaSet);
        Delay.SetValue(0.0);
        Delay.SetParameterSet(clParaSet);
    }
}

void CSPropPBCExcitation::SetPropagationDir(double val, int Component)
{
    if ((Component<0) || (Component>=3)) return;
    PropagationDir[Component].SetValue(val);
}

void CSPropPBCExcitation::SetPropagationDir(const std::string val, int Component)
{
    if ((Component<0) || (Component>=3)) return;
    PropagationDir[Component].SetValue(val);
}

double CSPropPBCExcitation::GetPropagationDir(int Component)
{
    if ((Component<0) || (Component>=3)) return 0;
    return PropagationDir[Component].GetValue();
}

const std::string CSPropPBCExcitation::GetPropagationDirString(int Comp)
{
    if ((Comp<0) || (Comp>=3)) return NULL;
    return PropagationDir[Comp].GetString();
}


bool CSPropPBCExcitation::Update(std::string *ErrStr, bool type)
{
    bool bOK=true;
    int EC=0;
    for (unsigned int i=0;i<3;++i)
    {
        if(type)
            EC = SINExcitation[i].Evaluate();
        if(!type)
            EC = COSExcitation[i].Evaluate();
        if (EC!=ParameterScalar::NO_ERROR) bOK=false;
        if ((EC!=ParameterScalar::NO_ERROR)  && (ErrStr!=NULL))
        {
            std::stringstream stream;
            stream << std::endl << "Error in PBC Excitation-Property Excitaion-Value (ID: " << uiID << "): ";
            ErrStr->append(stream.str());
            PSErrorCode2Msg(EC,ErrStr);
        }
        EC=PropagationDir[i].Evaluate();
        if (EC!=ParameterScalar::NO_ERROR) bOK=false;
        if ((EC!=ParameterScalar::NO_ERROR)  && (ErrStr!=NULL))
        {
            std::stringstream stream;
            stream << std::endl << "Error in PBC Excitation-Property PropagationDir-Value (ID: " << uiID << "): ";
            ErrStr->append(stream.str());
            PSErrorCode2Msg(EC,ErrStr);
        }
    }
    EC=m_Frequency.Evaluate();
    if (EC!=ParameterScalar::NO_ERROR) bOK=false;
    if ((EC!=ParameterScalar::NO_ERROR)  && (ErrStr!=NULL))
    {
        std::stringstream stream;
        stream << std::endl << "Error in PBC Excitation-Property Frequency-Value";
        ErrStr->append(stream.str());
        PSErrorCode2Msg(EC,ErrStr);
    }
    EC=Delay.Evaluate();
    if (EC!=ParameterScalar::NO_ERROR) bOK=false;
    if ((EC!=ParameterScalar::NO_ERROR)  && (ErrStr!=NULL))
    {
        std::stringstream stream;
        stream << std::endl << "Error in PBC Excitation-Property Delay-Value";
        ErrStr->append(stream.str());
        PSErrorCode2Msg(EC,ErrStr);
    }
    return bOK;
}

bool CSPropPBCExcitation::Write2XML(TiXmlNode& root, bool parameterised, bool sparse)
{
    if (CSProperties::Write2XML(root,parameterised,sparse) == false) return false;
    TiXmlElement* prop=root.ToElement();
    if (prop==NULL) return false;

    prop->SetAttribute("Number",(int)uiNumber);
    WriteTerm(m_Frequency,*prop,"Frequency",parameterised);
    WriteTerm(Delay,*prop,"Delay",parameterised);

    prop->SetAttribute("Type",iExcitType);
    WriteVectorTerm(SINExcitation,*prop,"SINExcite",parameterised);
    WriteVectorTerm(COSExcitation,*prop,"COSExcite",parameterised);
    TiXmlElement COSWeight("COSWeight");
    TiXmlElement SINWeight("SINWeight");
    WriteTerm(COSWeightFct[0],COSWeight,"X",parameterised);
    WriteTerm(COSWeightFct[1],COSWeight,"Y",parameterised);
    WriteTerm(COSWeightFct[2],COSWeight,"Z",parameterised);
    WriteTerm(SINWeightFct[0],SINWeight,"X",parameterised);
    WriteTerm(SINWeightFct[1],SINWeight,"Y",parameterised);
    WriteTerm(SINWeightFct[2],SINWeight,"Z",parameterised);
    prop->InsertEndChild(COSWeight);
    prop->InsertEndChild(SINWeight);

    WriteVectorTerm(PropagationDir,*prop,"PropDir",parameterised);

    return true;
}

bool CSPropPBCExcitation::ReadFromXML(TiXmlNode &root)
{
    if (CSProperties::ReadFromXML(root)==false) return false;

    TiXmlElement *prop = root.ToElement();
    if (prop==NULL) return false;

    int iHelp;
    if (prop->QueryIntAttribute("Number",&iHelp)!=TIXML_SUCCESS) uiNumber=0;
    else uiNumber=(unsigned int)iHelp;

    if (prop->QueryIntAttribute("Type",&iExcitType)!=TIXML_SUCCESS) return false;

    if (ReadVectorTerm(SINExcitation,*prop,"Excite",0.0)==false) return false;
    ReadTerm(m_Frequency,*prop,"Frequency");
    ReadTerm(Delay,*prop,"Delay");

    TiXmlElement *SINweight = prop->FirstChildElement("SINWeight");
    TiXmlElement *COSweight = prop->FirstChildElement("COSWeight");
    if (SINweight!=NULL)
    {
        ReadTerm(SINWeightFct[0],*SINweight,"X");
        ReadTerm(SINWeightFct[1],*SINweight,"Y");
        ReadTerm(SINWeightFct[2],*SINweight,"Z");
    }
    if (COSweight!=NULL)
    {
        ReadTerm(COSWeightFct[0],*COSweight,"X");
        ReadTerm(COSWeightFct[1],*COSweight,"Y");
        ReadTerm(COSWeightFct[2],*COSweight,"Z");
    }
    ReadVectorTerm(PropagationDir,*prop,"PropDir",0.0);

    return true;
}

void CSPropPBCExcitation::ShowPropertyStatus(std::ostream& stream)
{
    CSProperties::ShowPropertyStatus(stream);
    stream << " --- PBC Excitation Properties --- " << std::endl;
    stream << "  Type: " << iExcitType << std::endl;
    stream << "  Active directions: " << ActiveDir[0] << "," << ActiveDir[1] << "," << ActiveDir[2] << std::endl;
    stream << "  SINExcitation\t: " << SINExcitation[0].GetValueString() << ", "  << SINExcitation[1].GetValueString() << ", "  << SINExcitation[2].GetValueString()  << std::endl;
    stream << "  SINWeighting\t: " << SINWeightFct[0].GetValueString() << ", "  << SINWeightFct[1].GetValueString() << ", "  << SINWeightFct[2].GetValueString()  << std::endl;
    stream << "  SINExcitation\t: " << COSExcitation[0].GetValueString() << ", "  << COSExcitation[1].GetValueString() << ", "  << COSExcitation[2].GetValueString()  << std::endl;
    stream << "  COSWeighting\t: " << COSWeightFct[0].GetValueString() << ", "  << COSWeightFct[1].GetValueString() << ", "  << COSWeightFct[2].GetValueString()  << std::endl;
    stream << "  Propagation Dir: " << PropagationDir[0].GetValueString() << ", "  << PropagationDir[1].GetValueString() << ", "  << PropagationDir[2].GetValueString()  << std::endl;
    stream << "  Delay\t\t: " << Delay.GetValueString() << std::endl;
}
