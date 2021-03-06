/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLCommandLineModuleNode.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.3 $

=========================================================================auto=*/

#ifndef __vtkMRMLCommandLineModuleNode_h
#define __vtkMRMLCommandLineModuleNode_h

/// MRML includes
#include <vtkMRMLNode.h>

#include "vtkMRMLCLIWin32Header.h"

class ModuleDescription;

/// \brief MRML node for representing the parameters allowing to run a command line module
class VTK_MRML_CLI_EXPORT vtkMRMLCommandLineModuleNode : public vtkMRMLNode
{
public:
  static vtkMRMLCommandLineModuleNode *New();
  vtkTypeMacro(vtkMRMLCommandLineModuleNode, vtkMRMLNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual vtkMRMLNode* CreateNodeInstance();

  /// Set node attributes
  virtual void ReadXMLAttributes(const char** atts);

  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  /// Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node);

  /// Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName()
    {return "CommandLineModule";}

  /// Get/Set the module description object. THe module description
  /// object is used to cache the current settings for the module.
  const ModuleDescription& GetModuleDescription() const;
  ModuleDescription& GetModuleDescription();
  void SetModuleDescription(const ModuleDescription& description);

  typedef enum { Idle=0, Scheduled=1, Running=2, Completed=3, CompletedWithErrors=4, Cancelled=5 } StatusType;

  /// Set the status of the node (Idle, Scheduled, Running,
  /// Completed).  The "modify" parameter indicates whether the object
  /// can be modified by the call.
  void SetStatus(StatusType status, bool modify=true);
  StatusType GetStatus() const;
  const char* GetStatusString() const;

  /// Read a parameter file. This will set any parameters that
  /// parameters in this ModuleDescription.
  bool ReadParameterFile(const std::string& filename);

  /// Write a parameter file. This will output any parameters that
  /// parameters in this
  /// ModuleDescription. "withHandlesToBulkParameters" allows to
  /// control whether all parameters are written or just the
  /// parameters with simple IO mechanisms.
  bool WriteParameterFile(const std::string& filename, bool withHandlesToBulkParameters = true);

  /// Get/Set a parameter for the module.
  bool SetParameterAsString(const char* name, const std::string& value);
  bool SetParameterAsInt(const char* name, int value);
  bool SetParameterAsBool(const char* name, bool value);
  bool SetParameterAsDouble(const char* name, double value);
  bool SetParameterAsFloat(const char* name, float value);

  std::string GetParameterAsString(const char* name) const;

  void SetModuleDescription(const char *name);
  std::string GetModuleVersion() const;
  std::string GetModuleTitle() const;
  std::string GetModuleTarget() const;
  std::string GetModuleType() const;
  bool IsValidGroupId(unsigned int group) const;
  bool IsValidParamId(unsigned int group, unsigned int param) const;
  unsigned int GetNumberOfParameterGroups() const;
  unsigned int GetNumberOfParametersInGroup(unsigned int group) const;
  std::string GetParameterGroupLabel(unsigned int group) const;
  std::string GetParameterGroupDescription(unsigned int group) const;
  std::string GetParameterGroupAdvanced(unsigned int group) const;
  std::string GetParameterTag(unsigned int group, unsigned int param) const;
  std::string GetParameterType(unsigned int group, unsigned int param) const;
  std::string GetParameterArgType(unsigned int group, unsigned int param) const;
  std::string GetParameterName(unsigned int group, unsigned int param) const;
  std::string GetParameterLongFlag(unsigned int group, unsigned int param) const;
  std::string GetParameterLabel(unsigned int group, unsigned int param) const;
  std::string GetParameterConstraints(unsigned int group, unsigned int param) const;
  std::string GetParameterMaximum(unsigned int group, unsigned int param) const;
  std::string GetParameterMinimum(unsigned int group, unsigned int param) const;
  std::string GetParameterDescription(unsigned int group, unsigned int param) const;
  std::string GetParameterChannel(unsigned int group, unsigned int param) const;
  std::string GetParameterIndex(unsigned int group, unsigned int param) const;
  std::string GetParameterDefault(unsigned int group, unsigned int param) const;
  std::string GetParameterFlag(unsigned int group, unsigned int param) const;
  std::string GetParameterMultiple(unsigned int group, unsigned int param) const;
  std::string GetParameterFileExtensions(unsigned int group, unsigned int param) const;
  std::string GetParameterCoordinateSystem(unsigned int group, unsigned int param) const;

  /// Methods to manage the master list of module description prototypes
  static int GetNumberOfRegisteredModules();
  static const char* GetRegisteredModuleNameByIndex(int idx);
  static void RegisterModuleDescription(ModuleDescription md);
  static bool HasRegisteredModule(const std::string& name);
  static ModuleDescription GetRegisteredModuleDescription(const std::string& name);

protected:
  void AbortProcess();

private:
  vtkMRMLCommandLineModuleNode();
  ~vtkMRMLCommandLineModuleNode();
  vtkMRMLCommandLineModuleNode(const vtkMRMLCommandLineModuleNode&);
  void operator=(const vtkMRMLCommandLineModuleNode&);

  class vtkInternal;
  vtkInternal * Internal;
};

#endif
