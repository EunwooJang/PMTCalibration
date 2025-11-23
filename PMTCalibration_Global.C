R__LOAD_LIBRARY(build/libRawObjs.dylib)

char* PMTCalib_Find_Dir_Name(const char* keywords_str) {

  DIR* dir = opendir(Root_Directory.c_str());
  struct dirent* entry;
  char* result = NULL;
    
  while ((entry = readdir(dir)) != NULL) {
    if (entry->d_type == DT_DIR && strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0) {
      if (strstr(entry->d_name, keywords_str) != NULL) {
        char subdir_path[512];
        snprintf(subdir_path, sizeof(subdir_path), "%s/%s", Root_Directory.c_str(), entry->d_name);
                
        result = (char*)malloc(strlen(subdir_path) + 2);
        snprintf(result, strlen(subdir_path) + 2, "%s/", subdir_path);
        closedir(dir);
        return result;
      }
    }
  }
  closedir(dir);
  return result;
}


char* PMTCalib_Find_File_Name(const char* keywords_str, int Run_Id = 0) {

  DIR* dir = opendir(Root_Directory.c_str());
  struct dirent* entry;
  char* result = NULL;

  while ((entry = readdir(dir)) != NULL) {
    if (entry->d_type == DT_DIR && strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0) {
      if (strstr(entry->d_name, keywords_str) != NULL) {
        char subdir_path[512];
        snprintf(subdir_path, sizeof(subdir_path), "%s/%s", Root_Directory.c_str(), entry->d_name);

        DIR* subdir = opendir(subdir_path);
        struct dirent* subentry;

        while ((subentry = readdir(subdir)) != NULL) {
          if (subentry->d_type == DT_REG && strstr(subentry->d_name, "RUN_") == subentry->d_name) {
            char* dot = strchr(subentry->d_name, '.');
            if (dot != NULL) {
              char* last_dot = strrchr(subentry->d_name, '.');
              if (last_dot != dot) {
                char* last_part = last_dot + 1;
                if (strlen(last_part) == 5) {
                  int all_digits = 1;
                  for (int i = 0; i < 5; i++) {
                    if (last_part[i] < '0' || last_part[i] > '9') {
                      all_digits = 0;
                      break;
                    }
                  }
                  if (all_digits) {
                    int file_id = atoi(last_part);
                    if (file_id == Run_Id) {
                      result = (char*)malloc(strlen(subdir_path) + strlen(subentry->d_name) + 2);
                      snprintf(result, strlen(subdir_path) + strlen(subentry->d_name) + 2, "%s/%s", subdir_path, subentry->d_name);
                      closedir(subdir);
                      closedir(dir);
                      return result;
                    }
                  }
                }
              }
            }
          }
        }
        closedir(subdir);
      }
    }
  }
  closedir(dir);
  return result;
}

void PMTCalib_Delete_Canvas() {
  delete gROOT->FindObject("h");
  delete gROOT->FindObject("can");
}

void PMTCalib_Save_Image(const char* file_path, TCanvas* can) {
  can->SaveAs(file_path);
}
