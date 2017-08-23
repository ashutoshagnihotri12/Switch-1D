function [ AL, AR, BL, BR ] = StateValueEuler(Q, QL, QLL, QR, QRR, limiter)
    BL = Q + 0.5*limiter((Q-QL)./(QR-Q+eps)).*(QR-Q);          %j+1/2-left-state
    BR = QR - 0.5*limiter((QR-Q)./(QRR-QR+eps)).*(QRR-QR);     %j+1/2-right-state
    AL = QL + 0.5*limiter((QL-QLL)./(Q-QL+eps)).*(Q-QL);       %j-1/2-left-state
    AR = Q - 0.5*limiter((Q-QL)./(QR-Q+eps)).*(QR-Q);          %j-1/2-right-state
end

