import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Fusion

Item {
    id: root
    property alias text: leftBtn.text
    property alias checked: leftBtn.checked
    function toggle() {leftBtn.toggle();toggled()}
    signal toggled
    signal deleteClicked

    implicitWidth: leftBtn.width + rightBtn.width + leftBtn.anchors.leftMargin + rightBtn.anchors.rightMargin
    implicitHeight: leftBtn.height
    states: [
        State {
            name: "randomBarcode"
            PropertyChanges {target: leftBtn; palette.button: "blue"}
        },
        State {
            name: "barcodeList"
            PropertyChanges {target: leftBtn; palette.button: "green"}
        },
        State {
            name: "constantRegion"
            PropertyChanges {target: leftBtn; palette.button: "orange"}
        }
    ]

    Button {
        id: leftBtn
        checkable: !checked
        width: contentItem.implicitWidth + leftPadding + rightPadding
        height: fontMetrics.height + topPadding + bottomPadding
        palette {
            buttonText: "black"
        }
        z: 1
        onToggled: root.toggled()

        FontMetrics {
            id: fontMetrics
        }
    }
    Button {
        id: rightBtn
        anchors.left: leftBtn.right
        anchors.leftMargin: -3
        anchors.verticalCenter: leftBtn.verticalCenter
        height: leftBtn.height
        leftPadding: rightPadding + 3
        width: contentItem.implicitWidth + leftPadding + rightPadding
        text: "❌"
        palette: leftBtn.palette
        z: 0
        onClicked: root.deleteClicked()
    }
}
