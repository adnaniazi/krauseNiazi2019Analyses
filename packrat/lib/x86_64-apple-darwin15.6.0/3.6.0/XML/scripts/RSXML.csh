if(`test -n "-L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/lib -lxml2 -lz -lpthread -licucore -lm"`) then

if(${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:-L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/lib -lxml2 -lz -lpthread -licucore -lm
else
   setenv LD_LIBRARY_PATH -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/lib -lxml2 -lz -lpthread -licucore -lm
endif

endif
