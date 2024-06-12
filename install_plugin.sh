echo "******** CLEAN ********"
mvn clean

echo -e "\n******** PACKAGE ********"
mvn package

echo -e "\n******** INSTALL ********"
cp "./target/Process_Pixels-0.1.0-SNAPSHOT.jar" "/c/ImageJ/plugins"
