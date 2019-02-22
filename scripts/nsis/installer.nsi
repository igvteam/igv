Name "IGV"

OutFile "IGV_Win_@VERSION@-installer.exe"
InstallDir "$PROGRAMFILES64\IGV_@VERSION@"

ShowInstDetails nevershow
ShowUninstDetails nevershow
SetCompressor /solid lzma
AutoCloseWindow true
Icon "IGV_@VERSION@\IGV_64.ico"
LicenseData license_win.txt
LicenseForceSelection radiobuttons

Page license
Page directory
Page instfiles
UninstPage instfiles

section
     setOutPath "$INSTDIR"
     File /a /r IGV_@VERSION@\*.*
     createShortCut "$DESKTOP\IGV_@VERSION@.lnk" "$INSTDIR\igv.bat" "" "$INSTDIR\IGV_64.ico"
     createDirectory "$SMPROGRAMS\IGV_@VERSION@"
     createShortCut "$SMPROGRAMS\IGV_@VERSION@\IGV.lnk" "$INSTDIR\igv.bat" "" "$INSTDIR\IGV_64.ico"
     #createShortCut "$SMPROGRAMS\IGV_@VERSION@\IGVTools.lnk" "$INSTDIR\igvtools_gui.bat"
     
     WriteUninstaller $INSTDIR\uninstaller.exe
     createShortCut "$SMPROGRAMS\IGV_@VERSION@\uninstaller.lnk" "$INSTDIR\uninstaller.exe"
sectionEnd

Function un.onInit
    MessageBox MB_YESNO "This will uninstall IGV_@VERSION@.  Continue?" IDYES NoAbort
      Abort ; causes uninstaller to quit.
    NoAbort:
FunctionEnd

#RequestExecutionLevel admin

section "Uninstall"
	setAutoClose true
	RMDir /r "$SMPROGRAMS\IGV_@VERSION@"
	Delete "$Desktop\IGV_@VERSION@.lnk"
	
	# NSIS bset-practice recommends not using RMDir /r $INSTDIR... 
	RMDir /r /REBOOTOK $INSTDIR\*.*
	RMDir /REBOOTOK $INSTDIR
sectionEnd