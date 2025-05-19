import QtQuick
import QtQuick.Controls
import Qt.labs.platform as Platform

pragma NativeMethodBehavior: AcceptThisObject


Pane {
    padding: 0
    ListView {
        id: listView
        signal deleteClicked(index: int)
        objectName: "listView"
        anchors.fill: parent
        clip: true
        model: myModel
        delegate: Section {
            id: delegate
            text: display
            state: typename
            checked: ListView.isCurrentItem
            onToggled: {
                listView.currentIndex = index
            }
            onDeleteClicked: listView.deleteClicked(index)
            transform: Scale{
                    id: transScale
                    xScale: 1
            }
            property alias xScale: transScale.xScale
        }
        orientation: ListView.Horizontal
        spacing: 2
        boundsBehavior: Flickable.StopAtBounds
        ScrollBar.horizontal: ScrollBar{}
        property int _oldcount: 0
        onCountChanged: {
            if (count > _oldcount && count - 1 != currentIndex)
                itemAtIndex(count - 1).toggle();
            _oldcount = count
        }
        add: Transition {
            NumberAnimation {property: "opacity"; from: 0; to: 1; duration: 150}
            NumberAnimation {property: "xScale"; from: 0; to: 1; duration: 150}
        }
        remove: Transition {
            NumberAnimation {property: "opacity"; to: 0; duration: 150}
            NumberAnimation {property: "xScale"; to: 0; duration: 150}
        }
        removeDisplaced: Transition {
            NumberAnimation {property: "x"; duration: 150}
        }
    }
}
