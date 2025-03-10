//
// Created by Lei on 4/4/2024.
//

#include "FileDialog.hpp"

NAMESPACE_BEGIN(PCO)
namespace viewer {

    std::string fileDialogOpen() {
        const int FILE_DIALOG_MAX_BUFFER = 1024;
        char buffer[FILE_DIALOG_MAX_BUFFER];
        buffer[0] = '\0';
        buffer[FILE_DIALOG_MAX_BUFFER - 1] = 'x'; // Initialize last character with a char != '\0'

#ifdef __APPLE__
        // For apple use applescript hack
  FILE * output = popen(
    "osascript -e \""
    "   tell application \\\"System Events\\\"\n"
    "           activate\n"
    "           set existing_file to choose file\n"
    "   end tell\n"
    "   set existing_file_path to (POSIX path of (existing_file))\n"
    "\" 2>/dev/null | tr -d '\n' ","r");
  if (output)
  {
    auto ret = fgets(buffer, FILE_DIALOG_MAX_BUFFER, output);
    if (ret == NULL || ferror(output))
    {
      // I/O error
      buffer[0] = '\0';
    }
    if (buffer[FILE_DIALOG_MAX_BUFFER - 1] == '\0')
    {
      // File name too long, buffer has been filled, so we return empty string instead
      buffer[0] = '\0';
    }
  }
#elif defined _WIN32

        // Use native windows file dialog box
        // (code contributed by Tino Weinkauf)

        OPENFILENAME ofn;       // common dialog box structure
        char szFile[260];       // buffer for file name

        // Initialize OPENFILENAME
        ZeroMemory(&ofn, sizeof(ofn));
        ofn.lStructSize = sizeof(ofn);
        ofn.hwndOwner = NULL;
        ofn.lpstrFile = szFile;
        // Set lpstrFile[0] to '\0' so that GetOpenFileName does not
        // use the contents of szFile to initialize itself.
        ofn.lpstrFile[0] = '\0';
        ofn.nMaxFile = sizeof(szFile);
        ofn.lpstrFilter = "*.*\0";//off\0*.off\0obj\0*.obj\0mp\0*.mp\0";
        ofn.nFilterIndex = 1;
        ofn.lpstrFileTitle = NULL;
        ofn.nMaxFileTitle = 0;
        ofn.lpstrInitialDir = NULL;
        ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

        // Display the Open dialog box.
        int pos = 0;
        if (GetOpenFileName(&ofn) == TRUE) {
            while (ofn.lpstrFile[pos] != '\0') {
                buffer[pos] = (char) ofn.lpstrFile[pos];
                pos++;
            }
        }
        buffer[pos] = 0;
#else

        // For linux use zenity
  FILE * output = popen("/usr/bin/zenity --file-selection","r");
  if (output)
  {
    auto ret = fgets(buffer, FILE_DIALOG_MAX_BUFFER, output);
    if (ret == NULL || ferror(output))
    {
      // I/O error
      buffer[0] = '\0';
    }
    if (buffer[FILE_DIALOG_MAX_BUFFER - 1] == '\0')
    {
      // File name too long, buffer has been filled, so we return empty string instead
      buffer[0] = '\0';
    }
  }

  // Replace last '\n' by '\0'
  if(strlen(buffer) > 0)
  {
    buffer[strlen(buffer)-1] = '\0';
  }

#endif
        return std::string(buffer);
    }

    std::string fileDialogSave() {
        const int FILE_DIALOG_MAX_BUFFER = 1024;
        char buffer[FILE_DIALOG_MAX_BUFFER];
        buffer[0] = '\0';
        buffer[FILE_DIALOG_MAX_BUFFER - 1] = 'x'; // Initialize last character with a char != '\0'

#ifdef __APPLE__
        // For apple use applescript hack
  // There is currently a bug in Applescript that strips extensions off
  // of chosen existing files in the "choose file name" dialog
  // I'm assuming that will be fixed soon
  FILE * output = popen(
    "osascript -e \""
    "   tell application \\\"System Events\\\"\n"
    "           activate\n"
    "           set existing_file to choose file name\n"
    "   end tell\n"
    "   set existing_file_path to (POSIX path of (existing_file))\n"
    "\" 2>/dev/null | tr -d '\n' ","r");
  if (output)
  {
    auto ret = fgets(buffer, FILE_DIALOG_MAX_BUFFER, output);
    if (ret == NULL || ferror(output))
    {
      // I/O error
      buffer[0] = '\0';
    }
    if (buffer[FILE_DIALOG_MAX_BUFFER - 1] == '\0')
    {
      // File name too long, buffer has been filled, so we return empty string instead
      buffer[0] = '\0';
    }
  }
#elif defined _WIN32

        // Use native windows file dialog box
        // (code contributed by Tino Weinkauf)

        OPENFILENAME ofn;       // common dialog box structure
        char szFile[260];       // buffer for file name

        // Initialize OPENFILENAME
        ZeroMemory(&ofn, sizeof(ofn));
        ofn.lStructSize = sizeof(ofn);
        ofn.hwndOwner = NULL;//hwnd;
        ofn.lpstrFile = szFile;
        // Set lpstrFile[0] to '\0' so that GetOpenFileName does not
        // use the contents of szFile to initialize itself.
        ofn.lpstrFile[0] = '\0';
        ofn.nMaxFile = sizeof(szFile);
        ofn.lpstrFilter = "";
        ofn.nFilterIndex = 1;
        ofn.lpstrFileTitle = NULL;
        ofn.nMaxFileTitle = 0;
        ofn.lpstrInitialDir = NULL;
        ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

        // Display the Open dialog box.
        int pos = 0;
        if (GetSaveFileName(&ofn) == TRUE) {
            while (ofn.lpstrFile[pos] != '\0') {
                buffer[pos] = (char) ofn.lpstrFile[pos];
                pos++;
            }
            buffer[pos] = 0;
        }

#else
        // For every other machine type use zenity
  FILE * output = popen("/usr/bin/zenity --file-selection --save","r");
  if (output)
  {
    auto ret = fgets(buffer, FILE_DIALOG_MAX_BUFFER, output);
    if (ret == NULL || ferror(output))
    {
      // I/O error
      buffer[0] = '\0';
    }
    if (buffer[FILE_DIALOG_MAX_BUFFER - 1] == '\0')
    {
      // File name too long, buffer has been filled, so we return empty string instead
      buffer[0] = '\0';
    }
  }

  // Replace last '\n' by '\0'
  if(strlen(buffer) > 0)
  {
    buffer[strlen(buffer)-1] = '\0';
  }

#endif
        return std::string(buffer);
    }

} // namespace viewer
NAMESPACE_END(PCO)