import ROOT
import sys

def check_workspace(root_file_path, ws_name="wsMy", pdf_name="massModel", dataset_name="ds_mass"):
    """
    ROOT 파일 안의 작업공간 내용물을 검증하는 스크립트.
    """
    print(f"--- 검증 시작: {root_file_path} ---")
    
    tfile = ROOT.TFile.Open(root_file_path)
    if not tfile or tfile.IsZombie():
        print("결과: [실패] ROOT 파일을 열 수 없습니다.")
        return

    workspace = tfile.Get(ws_name)
    if not workspace:
        print(f"결과: [실패] 파일 안에 '{ws_name}' 이름의 작업공간(RooWorkspace)을 찾을 수 없습니다.")
        tfile.Close()
        return

    print(f"'{ws_name}' 작업공간을 성공적으로 불러왔습니다. 내용물을 확인합니다...")
    
    pdf = workspace.pdf(pdf_name)
    data = workspace.data(dataset_name)
    nsig = workspace.var("nSig")
    nbkg = workspace.var("nBkg")

    all_found = True
    if not pdf:
        print(f" - PDF '{pdf_name}': [실패]")
        all_found = False
    else:
        print(f" - PDF '{pdf_name}': [성공]")

    if not data:
        print(f" - 데이터셋 '{dataset_name}': [실패]")
        all_found = False
    else:
        print(f" - 데이터셋 '{dataset_name}': [성공]")

    if not nsig or not nbkg:
        print(f" - 수율 변수 'nSig'/'nBkg': [실패]")
        all_found = False
    else:
        print(f" - 수율 변수 'nSig'/'nBkg': [성공]")

    print("--- 검증 완료 ---")
    if all_found:
        print("최종 결과: [성공] 작전 기반이 완벽하게 확보되었습니다.")
    else:
        print("최종 결과: [실패] sPlot에 필요한 일부 요소가 누락되었습니다.")

    tfile.Close()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        check_workspace(sys.argv[1])
    else:
        print("사용법: python check_workspace.py [검사할 .root 파일 경로]")