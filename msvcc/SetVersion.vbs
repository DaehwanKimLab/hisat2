
' ask the scripting runtime environemt for access to files
set FS = CreateObject("Scripting.FileSystemObject")
set FileHandler = FS.GetFile("../VERSION")
set inputTextStream = FileHandler.OpenAsTextStream(1)
version = inputTextStream.ReadLine

WScript.Echo  "Version:  " & version

set outputTextStream = FS.CreateTextFile("../version.h",true)
outputTextStream.WriteLine "#define HISAT2_VERSION """ & version & """"
outputTextStream.close


